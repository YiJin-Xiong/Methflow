# Set the base image to Ubuntu and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags
FROM docker.1ms.run/nvidia/cuda:12.3.2-cudnn9-devel-ubuntu20.04

# Author and maintainer
LABEL description="Methflow" \
      author="Yijin Xiong <yshico@163.com>"

ARG DNAME="Methflow"

# Add CUDA to PATH
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=$PATH:$CUDA_HOME/bin

ARG BUILD_PACKAGES="wget apt-transport-https procps git curl git-lfs build-essential"
ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -q && \
    DEBIAN_FRONTEND=${DEBIAN_FRONTEND} && \
    apt-get install -q --yes ${BUILD_PACKAGES} && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#Install miniconda
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh && \
    /bin/bash Miniconda.sh -b -p /opt/conda && \
    rm Miniconda.sh

# Adding conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Install dorado
ARG DORADO_VERSION="0.9.5"
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz && \
    tar -xzf dorado-${DORADO_VERSION}-linux-x64.tar.gz \
    && mv dorado-${DORADO_VERSION}-linux-x64 /opt/dorado \
    && rm dorado-${DORADO_VERSION}-linux-x64.tar.gz

# Add dorado to PATH
ENV PATH /opt/dorado/bin:$PATH
RUN dorado --version

# Install rockfish
RUN git clone -b r10.4.1 https://github.com/lbcb-sci/rockfish.git --single-branch /opt/rockfish && \
    cd /opt/rockfish && \
    pip install --extra-index-url https://download.pytorch.org/whl/cu123 .
RUN rockfish

# Install deepsignal3
RUN git clone https://github.com/PengNi/deepsignal3.git /opt/deepsignal3 && \
    cd /opt/deepsignal3 && \
    python setup.py install
RUN deepsignal3 --version

# Create the environment
COPY environment.yml /
RUN conda env create --name ${DNAME} --file=environment.yml && conda clean -a

## install clair3 environment
## https://github.com/HKU-BAL/Clair3/blob/main/Dockerfile
#COPY environment-clair3.yml /
#RUN conda env create --name clair3 --file=environment-clair3.yml && conda clean -a

# Make RUN commands use the new environment
# name need to be the same with the above ${DNAME}
SHELL ["conda", "run", "-n", "Methflow", "/bin/bash", "-c"]

# clear pip cache
RUN pip cache purge

# download ccsmeth model
# RUN echo '199.232.28.133 raw.githubusercontent.com' | tee -a /etc/hosts && \
#     mkdir -p /opt/models/ccsmeth && \
#     cd /opt/models/ccsmeth && \
#     wget -q https://github.com/PengNi/basemods-models/raw/master/ccsmeth/model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v2.ckpt && \
#     wget -q https://github.com/PengNi/basemods-models/raw/master/ccsmeth/model_ccsmeth_5mCpG_aggregate_attbigru_b11.v2p.ckpt && \
#     ls -lh

# Set env path into PATH
ENV PATH /opt/conda/envs/${DNAME}/bin:$PATH
USER root
WORKDIR /data/

RUN cd /data

CMD ["bash"]

