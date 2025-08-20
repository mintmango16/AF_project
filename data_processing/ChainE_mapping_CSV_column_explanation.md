**uniprot\_pos**

UniProt Q6PSU2의 1-based 잔기 번호(1..172).

이 값이 있고 pdb\_\*가 None이면 구조가 없는 자리(모델링되지 않은 잔기) 를 뜻합니다.



**uniprot\_aa**

해당 UniProt 위치의 1-letter 아미노산.

IEDB 에피토프 레이블 등 서열 기반 레이블을 여기에 붙이면 됩니다.



**pdb\_chain**

매핑된 PDB의 체인 ID (이 경우, chain "E").

UniProt에는 없고 PDB에만 있는 행일 경우(드물게) uniprot\_\*가 None이고 여기가 채워집니다.



**pdb\_modeled\_index**

PDB 체인 내에서 “좌표가 있는 폴리펩타이드 잔기”만 세었을 때의 연속 인덱스(1..N).

즉, 결손/누락 잔기는 건너뛰고 실제로 구조가 있는 잔기만 1,2,3…으로 번호를 매긴 값입니다.

구조 피처(DSSP, RSA 등) 배열의 인덱스를 맞출 때 유용합니다.



**pdb\_resseq**

PDB 파일의 author residue number.

비연속/점프/음수/큰 번호 등 원본 PDB 표기 그대로라서, 사람이 보는 레퍼런스나 외부 도구 호환에 필요합니다.



**pdb\_icode**

insertion code(삽입코드). 없으면 빈 문자열.

PDB에서는 같은 번호에 A/B/...를 덧붙여 잔기를 삽입할 수 있어, resseq + icode가 유일 식별자가 됩니다.



**pdb\_aa**

PDB 잔기의 1-letter 아미노산(비표준은 보정: 예 MSE→M 등).

uniprot\_aa와 비교해 치환/변이/오타를 빠르게 확인할 수 있습니다.



**match**

uniprot\_aa == pdb\_aa이면 1, 다르면 0.

한쪽이 비어 있는(갭) 행은 None 입니다.

0이 발생하는 이유: 아이소폼/품종 차이, 비표준 보정, 입력 서열과의 미세 차이 등.



**행 해석 예시**

uniprot\_pos=85, uniprot\_aa='D', pdb\_modeled\_index=60, pdb\_resseq=120, pdb\_aa='D', match=1

→ UniProt 85번 D가 PDB 체인 E의 resseq 120(체인 내 60번째 모델링 잔기)와 일치.



uniprot\_pos=25, uniprot\_aa='Q', pdb\_\* = None

→ UniProt 25번은 8DB4에 좌표가 없음(누락/삭제 구간). 구조 피처는 NaN/마스크 처리.



uniprot\_pos=None, pdb\_modeled\_index=5, pdb\_resseq=42, pdb\_aa='S'

→ PDB에만 존재하는 자리(정렬상 갭). 대개 드물고, 보통은 데이터 병합에서 제외합니다.



match=0

→ 서열 불일치. 실제 변이/비표준 보정/서열 출처 차이를 의심하고 확인.

