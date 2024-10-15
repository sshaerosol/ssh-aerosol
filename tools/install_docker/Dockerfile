FROM ubuntu:24.04

RUN apt-get -y update && apt-get install -y \
    python3 \
    python3-pip \
    cmake \
    vim \
    wget \
    unzip \
    curl \
    build-essential \
    zlib1g-dev \
    libcurl4-openssl-dev \
    m4 \
    gfortran \
    scons \
    libnetcdff-dev

WORKDIR /root

COPY setup_ssh-aerosol.sh /root/setup_ssh-aerosol.sh

RUN chmod +x /root/setup_ssh-aerosol.sh

RUN /root/setup_ssh-aerosol.sh

COPY ssh-aerosol /root/ssh-aerosol
