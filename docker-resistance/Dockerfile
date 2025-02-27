# Usar la imagen base de Ubuntu 20.04
FROM ubuntu:20.04

# Establecer variables de entorno para evitar preguntas interactivas durante la instalación
ENV DEBIAN_FRONTEND=noninteractive

# Actualizar el sistema e instalar dependencias básicas
RUN apt-get update && apt-get install -y \
    wget \
    git \
    python3 \
    python3-pip \
    build-essential \
    libz-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libxslt1-dev \
    libffi-dev \
    libopenblas-dev \
    liblapack-dev \
    gfortran \
    libhdf5-serial-dev \
    libfreetype6-dev \
    libpng-dev \
    libopenmpi-dev \
    openmpi-bin \
    openmpi-common \
    openmpi-doc \
    && rm -rf /var/lib/apt/lists/*

# Instalar dependencias comunes de Python
RUN pip3 install --upgrade pip
RUN pip3 install biopython pandas numpy scipy matplotlib

# Instalar resfinder
RUN git clone https://bitbucket.org/genomicepidemiology/resfinder.git /opt/resfinder
RUN pip3 install cgecore tabulate

# Instalar virulencefinder
RUN git clone https://bitbucket.org/genomicepidemiology/virulencefinder.git /opt/virulencefinder
RUN pip3 install cgecore tabulate

# Instalar RGI (CARD Resistance Gene Identifier)
RUN git clone https://github.com/arpcard/rgi.git /opt/rgi
RUN pip3 install -r /opt/rgi/requirements.txt
RUN cd /opt/rgi && python3 setup.py install

# Instalar sccmec
RUN git clone https://github.com/egonozer/sccmec.git /opt/sccmec
RUN pip3 install biopython

# Establecer el directorio de trabajo
WORKDIR /workspace

# Comando por defecto al iniciar el contenedor
CMD ["/bin/bash"]
