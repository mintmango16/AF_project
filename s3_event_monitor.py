#!/usr/bin/env python3
"""
S3 Event Notifications + SQS를 이용한 실시간 파이프라인 자동화
"""

import boto3
import json
import logging
import subprocess
import time
import os
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pipeline_events.log'),
        logging.StreamHandler()
    ]
)

class S3EventPipelineMonitor:
    def __init__(self, config_file='config.json'):
        self.load_config(config_file)
        self.s3_client = boto3.client('s3')
        self.sqs_client = boto3.client('sqs')
        
    def load_config(self, config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        self.bucket_name = config['s3_bucket']
        self.sqs_queue_url = config['sqs_queue_url']
        self.nextflow_script = config['nextflow_script']
        self.work_dir = Path(config['work_directory'])
        self.output_prefix = config['output_prefix']
        
        self.work_dir.mkdir(exist_ok=True)
        logging.info(f"작업 디렉토리: {self.work_dir}")
    
    def process_sqs_messages(self):
        """SQS에서 S3 이벤트 메시지 처리"""
        try:
            response = self.sqs_client.receive_message(
                QueueUrl=self.sqs_queue_url,
                MaxNumberOfMessages=10,
                WaitTimeSeconds=20,  # Long polling
                MessageAttributeNames=['All']
            )
            
            messages = response.get('Messages', [])
            
            for message in messages:
                try:
                    # S3 이벤트 파싱
                    body = json.loads(message['Body'])
                    
                    # SNS를 통해 전달된 경우
                    if 'Message' in body:
                        s3_event = json.loads(body['Message'])
                    else:
                        s3_event = body
                    
                    # S3 레코드 처리
                    for record in s3_event.get('Records', []):
                        if record['eventName'].startswith('ObjectCreated'):
                            self.handle_s3_object_created(record)
                    
                    # 메시지 삭제 (처리 완료)
                    self.sqs_client.delete_message(
                        QueueUrl=self.sqs_queue_url,
                        ReceiptHandle=message['ReceiptHandle']
                    )
                    
                except Exception as e:
                    logging.error(f"메시지 처리 오류: {e}")
            
            return len(messages)
            
        except Exception as e:
            logging.error(f"SQS 메시지 수신 오류: {e}")
            return 0
    
    def handle_s3_object_created(self, record):
        """S3 객체 생성 이벤트 처리"""
        bucket = record['s3']['bucket']['name']
        key = record['s3']['object']['key']
        
        logging.info(f"S3 이벤트 수신: {bucket}/{key}")
        
        # FASTA 파일인지 확인
        if not (key.endswith('.fasta') or key.endswith('.fa')):
            logging.info(f"FASTA 파일이 아님, 무시: {key}")
            return
        
        if not key.startswith('input/'):
            logging.info(f"input/ 경로가 아님, 무시: {key}")
            return
        
        logging.info(f"새 FASTA 파일 감지: {key}")
        
        # 파이프라인 실행
        self.process_fasta_file(bucket, key)
    
    def process_fasta_file(self, bucket, s3_key):
        """FASTA 파일 처리 파이프라인"""
        file_id = Path(s3_key).stem
        timestamp = int(time.time())
        
        allergen_name = file_id.replace('_', ' ').replace('-', ' ').title()
        
        try:
            # 1. 파일 다운로드
            local_fasta = self.work_dir / f"{file_id}_{timestamp}.fasta"
            self.s3_client.download_file(bucket, s3_key, str(local_fasta))
            logging.info(f"파일 다운로드: {s3_key} -> {local_fasta}")
            
            # 2. 출력 디렉토리 설정
            output_dir = self.work_dir / f"output_{file_id}_{timestamp}"
            output_dir.mkdir(exist_ok=True)
            
            # 3. Nextflow 실행
            success = self.run_nextflow_pipeline(local_fasta, output_dir)
            
            if success:
                # 4. 결과 업로드
                self.upload_results(output_dir, f"{self.output_prefix}/{file_id}_{timestamp}")
                logging.info(f"파이프라인 완료: {file_id}")
                
                # 5. 메타데이터 생성
                self.create_metadata(file_id, s3_key, timestamp, True, allergen_name=allergen_name)
            else:
                logging.error(f"파이프라인 실패: {file_id}")
                self.create_metadata(file_id, s3_key, timestamp, False)
            
        except Exception as e:
            logging.error(f"파일 처리 오류: {e}")
            self.create_metadata(file_id, s3_key, timestamp, False, str(e))
    
    def run_nextflow_pipeline(self, fasta_file, output_dir):
        """Nextflow 파이프라인 실행"""
        try:
            cmd = [
                'nextflow', 'run', self.nextflow_script,
                '--fasta', str(fasta_file),
                '--outdir', str(output_dir),
                '--use_alphafold', 'true',
                '-resume'
            ]
            
            logging.info(f"Nextflow 시작: {fasta_file.name}")
            start_time = datetime.now()
            
            # 환경변수 설정
            env = os.environ.copy()
            env['NXF_OPTS'] = '-Xmx8g'  # Nextflow 메모리 제한
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200,  # 2시간
                env=env,
                cwd=self.work_dir
            )
            
            duration = datetime.now() - start_time
            
            if result.returncode == 0:
                logging.info(f"Nextflow 완료: {duration}")
                return True
            else:
                logging.error(f"Nextflow 실패: {result.stderr}")
                # 에러 로그를 파일로 저장
                with open(self.work_dir / 'nextflow_error.log', 'w') as f:
                    f.write(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            logging.error("Nextflow 타임아웃 (2시간)")
            return False
        except Exception as e:
            logging.error(f"Nextflow 실행 오류: {e}")
            return False
    
    def upload_results(self, local_dir, s3_prefix):
        """결과를 S3에 업로드"""
        try:
            uploaded = 0
            outputs_dir = local_dir / 'outputs'
            if outputs_dir.exists() and outputs_dir.is_dir():
                # outputs 폴더 내부 파일들을 직접 업로드
                upload_dir = outputs_dir
            else:
                # outputs 폴더가 없으면 전체 디렉토리 업로드
                upload_dir = local_dir
            
            for file_path in upload_dir.rglob('*'):
                if file_path.is_file():
                    # outputs 폴더를 건너뛰고 상대 경로 생성
                    if outputs_dir.exists():
                        relative_path = file_path.relative_to(outputs_dir)
                    else:
                        relative_path = file_path.relative_to(local_dir)
                        
                    s3_key = f"{s3_prefix}/{relative_path}"
                    
                    self.s3_client.upload_file(
                        str(file_path),
                        self.bucket_name,
                        s3_key
                    )
                    uploaded += 1
            
            logging.info(f"결과 업로드 완료: {uploaded}개 파일")
            
        except Exception as e:
            logging.error(f"업로드 오류: {e}")
    
    def create_metadata(self, file_id, s3_key, timestamp, success, error_msg=None, allergen_name=None):
        from datetime import datetime, timezone, timedelta
        
        # 한국 시간으로 변환
        kst = timezone(timedelta(hours=9))
        created_time = datetime.now(kst).strftime("%Y-%m-%dT%H:%M:%S+09:00")
        
        metadata = {
            'file_id': file_id,
            'allergen_name': allergen_name or file_id,
            'input_file': s3_key,
            'processing_time': created_time,  # 한국 시간 사용
            'timestamp': timestamp,
            'status': 'completed' if success else 'failed',
            'error': error_msg if error_msg else None
        }
        
        try:
            metadata_key = f"{self.output_prefix}/{file_id}_{timestamp}/metadata.json"
            self.s3_client.put_object(
                Bucket=self.bucket_name,
                Key=metadata_key,
                Body=json.dumps(metadata, indent=2)
            )
            logging.info(f"메타데이터 업로드: {metadata_key}")
            
        except Exception as e:
            logging.error(f"메타데이터 업로드 오류: {e}")
    
    def start_event_monitoring(self):
        """이벤트 기반 모니터링 시작"""
        logging.info("=== S3 Event 모니터링 시작 ===")
        logging.info(f"버킷: {self.bucket_name}")
        logging.info(f"SQS Queue: {self.sqs_queue_url}")
        logging.info(f"Nextflow 스크립트: {self.nextflow_script}")
        logging.info("=" * 50)
        
        while True:
            try:
                # SQS에서 메시지 수신 (Long Polling 20초)
                message_count = self.process_sqs_messages()
                
                if message_count == 0:
                    logging.debug("새 메시지 없음")
                else:
                    logging.info(f"{message_count}개 메시지 처리됨")
                
            except KeyboardInterrupt:
                logging.info("모니터링 중단")
                break
            except Exception as e:
                logging.error(f"모니터링 오류: {e}")
                time.sleep(10)  # 오류 시 10초 대기

if __name__ == "__main__":
    if not Path('config.json').exists():
        print("config.json 파일이 필요합니다.")
        exit(1)
    
    monitor = S3EventPipelineMonitor()
    monitor.start_event_monitoring()