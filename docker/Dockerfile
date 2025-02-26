FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Taller resistencia-INMEGEN"
LABEL version="1.0"

# Build the image: docker build -t laugoro/resistance-workshop-inmegen:public -f Dockerfile .
# To run the container: docker run -i -t laugoro/resistance-workshop-inmegen:public /bin/bash

# Actualiza el sistema operativo e instala las herramientas adecuadas
RUN apt-get update && \
    apt-get install -y \
        ant \
        automake \
        apt-utils \
        build-essential \
        bzip2 \
        ca-certificates \
        cutadapt \
        curl \
        gawk \
        gcc \
        git \
        g++ \
        hmmer \
        openjdk-8-jdk \
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
        python3-dev \
        python3-numpy \
        python3-setuptools \
        python3-pip \
        rna-star \
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

#### Descargar e instalar FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod a+x FastQC/fastqc && \
    ln -s /usr/src/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

#### Descargar e instalar MultiQC
RUN pip3 install --upgrade numpy
RUN pip3 install multiqc 



#### Descargar e instalar BWA
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xjvf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    ln -s /usr/src/bwa-0.7.17/bwa /usr/local/bin/bwa && \
    cd .. && \
    rm bwa-0.7.17.tar.bz2

#### Descargar e instalar Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar -xjvf samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm samtools-1.12.tar.bz2


#### Descargar e instalar mosdepth
RUN wget "https://github.com/brentp/mosdepth/releases/download/v0.3.6/mosdepth" && \
    cp ./mosdepth /usr/local/bin/ && \
    chmod a+x /usr/local/bin/mosdepth && \
    ln -s /usr/local/bin/mosdepth /usr/bin/mosdepth

#### Decargar e instalar Fastp
RUN wget http://opengene.org/fastp/fastp && \
    cp ./fastp /usr/local/bin/ && \
    chmod a+x /usr/local/bin/fastp && \
    ln -s /usr/local/bin/fastp /usr/bin/fastp

#### Descargar e instalar SPADES
RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz && \
    tar -xzf SPAdes-4.0.0-Linux.tar.gz && \
    rm SPAdes-4.0.0-Linux.tar.gz

#### Descargar e instalar QUAST
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.3.0.tar.gz && \
    tar -xzf quast-5.3.0.tar.gz && \
    rm quast-5.3.0.tar.gz 

ENV PATH="${PATH}:/usr/src/SPAdes-4.0.0-Linux/bin:/usr/src/quast-5.3.0"

# Ejecutar 
CMD ["/bin/bash"]

