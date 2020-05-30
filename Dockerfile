################## BASE IMAGE ######################

FROM ubuntu:16.04

################## METADATA ######################
LABEL base.image="ubuntu:16.04"
LABEL software="msg"
LABEL about.home="https://github.com/YourePrettyGood/msg"
LABEL about.documentation="https://github.com/YourePrettyGood/msg"
LABEL license="https://github.com/YourePrettyGood/msg/blob/master/LICENSE"
################## MAINTAINER ######################

MAINTAINER Yuan(Alvin) Chen <ychenbioinfo@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python          \
	    python-dev      \
        python-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        gfortran        \
        mercurial       \
        subversion      \
        libncurses5-dev \
        libncursesw5-dev \
        libbz2-dev      \
        liblzma-dev     \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN mkdir /msg /data /usr/opt

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/msg

RUN conda config --add channels r
RUN conda config --add channels bioconda

RUN conda upgrade conda

# install bwa - failed to compile 0.5.7 provided in the dependencies, use 0.5.9 from conda 
RUN conda install bwa=0.5.9-0 r r-r.methodss3 

RUN pip install --upgrade pip
RUN pip install cython==0.29.19 
RUN pip install pysam==0.15.4 biopython==1.76 numpy==1.16.6 mailer==0.8.1

COPY ./ /msg
COPY ./cmdline /msg/cmdline
COPY ./dependencies/zlib-1.2.5.tar.bz2 /usr/opt
COPY ./dependencies/samtools-0.1.9.tar.bz2 /usr/opt
COPY ./dependencies/Pyrex-0.9.9.tar.gz /usr/opt

WORKDIR /usr/opt

# install libz
RUN tar -jvxf zlib-1.2.5.tar.bz2 && \
    cd zlib-1.2.5 && \
    ./configure && \
    make && \
    make install

# install samtools
RUN tar -jvxf samtools-0.1.9.tar.bz2 && \ 
    cd samtools-0.1.9 && \
    make

# install Pyrex
RUN tar -zvxf Pyrex-0.9.9.tar.gz && \
    cd Pyrex-0.9.9 && \
    python setup.py install

ENV PATH /usr/opt/samtools-0.1.9:$PATH:

WORKDIR /msg

RUN gcc -O2 -lm -o countalleles countalleles.c && \
    gcc -O2 -lm -o hmmprobs hmmprobs.c 

RUN sh test_dependencies_R.sh

WORKDIR /data