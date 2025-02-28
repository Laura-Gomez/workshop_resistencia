FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Servicios INMEGEN"
LABEL version="1.0"

# Building takes several hours! Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t pipelinesinmegen/bacterias_inmegen:c1 -f Dockerfile .
# To run the container: docker run -i -t  pipelinesinmegen/bacterias_inmegen:c1 /bin/bash

# Actualiza el sistema operativo e instala las herramientas adecuadas
RUN apt-get update && \
    apt-get install -y \
        ant \
        automake \
        apt-utils \
        build-essential \
        bzip2 \
        ca-certificates \
        gawk \
        gcc \
        git \
        g++ \
        hmmer \
        openjdk-11-jdk \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libpng-dev \
        libpcre2-dev \ 
        libpcre3 \
        libpcre3-dev \
	libncurses5-dev \
        libnss-sss \	
        libcairo2-dev \
        libxt-dev \
        libfontconfig1-dev \
        libudunits2-dev \
        libgdal-dev \
        less \
        locales \
        make \
        nano \
	nim \
        parallel \
        perl \
        pigz \
        pkg-config \
        software-properties-common \
        tabix \
        tar \
        tzdata \
        unzip \
        vim \
        wget \
        zip \
        zlib1g \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src

#### Descargar e instalar Miniconda 3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Agregar la ruta de conda al PATH
ENV PATH /conda/bin:$PATH

# Crear y activar un ambiente conda por defecto
RUN conda create -y -n bact_env python=3.12
RUN echo "source activate bact_env" >> ~/.bashrc 
#RUN echo "pip install cgelib" >> ~/.bashrc
#RUN echo "pip install resfinder" >> ~/.bashrc

# Canales de conda
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

#RUN ./.bashrc

## instalar programas
RUN conda install conda-forge::biopython=1.84
RUN conda install conda-forge::tabulate=0.9.0
#RUN conda install --use-local --update-deps bioconda::resfinder
#RUN conda install bioconda::resfinder 
RUN pip install cgelib && pip install resfinder

RUN conda init && bash ~/.bashrc && . ~/.bashrc
RUN conda init && . /root/.bashrc && conda activate bact_env && pip install cgelib && pip install resfinder


RUN echo "source activate bact_env" >> ~/.bashrc 


## fix resfinder python /conda/lib/python3.12/site-packages/resfinder/run_resfinder.py
RUN mkdir -p /conda/envs/bact_env/bin
RUN echo '#!/usr/bin/env bash' > /conda/envs/bact_env/bin/resfinder && \
    echo 'python /conda/lib/python3.12/site-packages/resfinder/run_resfinder.py ${@}' >> /conda/envs/bact_env/bin/resfinder && \
    chmod a+x /conda/envs/bact_env/bin/resfinder && \
    cp /conda/envs/bact_env/bin/resfinder /usr/local/bin/ && \
    chmod a+x /usr/local/bin/resfinder && \
    ln -s /usr/local/bin/resfinder /usr/bin/resfinder 


## instalar blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xzf  ncbi-blast-2.16.0+-x64-linux.tar.gz

## instalar kma
RUN git clone https://bitbucket.org/genomicepidemiology/kma.git && \
    cd kma && \
    make

## instalar virfinder
RUN conda install bioconda::virulencefinder

## instalar rgi
#RUN conda install conda-forge::mamba
#RUN mamba install --channel conda-forge --channel bioconda --channel defaults rgi=5.1.1
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    cp diamond /usr/local/bin/ && \
    ln -s /usr/local/bin/diamond /usr/bin/diamond

RUN pip install git+https://github.com/arpcard/rgi.git

## instalar sccmec
RUN conda install bioconda::sccmec

## instalar prokka
RUN conda install bioconda::prokka

## instalar hmmscan
RUN conda install bioconda::hmmer

### Otras variables de entorno ###
ENV TZ "America/Mexico_City"

ENV PATH=$PATH:/usr/src/ncbi-blast-2.16.0+/bin:/usr/src/kma

# Ejecutar 
CMD ["/bin/bash"]
