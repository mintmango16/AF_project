FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-devel
ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /app

# 컴파일러/빌드툴 (openfold 확장 빌드용)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ninja-build git curl wget unzip ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# nvcc 위치 및 아키텍처(※ CUDA 11.6이 지원하는 범위로 고정: 8.6/8.0/7.5/7.0/6.1/5.2)
ENV CUDA_HOME=/usr/local/cuda
ENV TORCH_CUDA_ARCH_LIST="8.6;8.0;7.5;7.0;6.1;5.2"
ENV MAX_JOBS=4


# mkdssp (micromamba + bioconda). "mkdssp" 바이너리만 노출하고, env의 Python은 PATH에 추가하지 않음.
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN set -eux; \
    curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
      | tar -xvj -C /usr/local -f - bin/micromamba; \
    /usr/local/bin/micromamba create -y -r ${MAMBA_ROOT_PREFIX} -n bio -c conda-forge -c bioconda dssp; \
    ln -s ${MAMBA_ROOT_PREFIX}/envs/bio/bin/mkdssp /usr/local/bin/mkdssp

# 파이썬 의존성은 반드시 베이스 파이썬(/opt/conda)로 설치
COPY requirements.txt /app/requirements.txt
RUN /opt/conda/bin/python -m pip install --upgrade pip && \
    /opt/conda/bin/python -m pip install --no-build-isolation --no-cache-dir -r requirements.txt

# GraphBepi 학습된 모델 가중치 복사
COPY model/my-foodepitope/EpiDot.ckpt /opt/models/my-foodepitope/


COPY . /app
