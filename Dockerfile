################## BASE IMAGE ######################

FROM ubuntu:18.04

################## METADATA ######################
LABEL base.image="ubuntu:18.04"
LABEL software="msg"
LABEL about.home="https://github.com/YourePrettyGood/msg"
LABEL about.documentation="https://github.com/YourePrettyGood/msg"
LABEL license="https://github.com/YourePrettyGood/msg/blob/master/LICENSE"
################## MAINTAINER ######################

MAINTAINER Yuan(Alvin) Chen <ychenbioinfo@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt bionic main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt bionic-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt bionic-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt bionic-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
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
        zlib1g-dev      \
        python-pyrex    \
        libswitch-perl  \
        ghostscript  && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.7.12.1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN mkdir /msg /data /usr/opt

ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/msg

RUN pip install --upgrade pip
#Addressing a conda bug:
RUN pip install tqdm
#The rest are MSG dependencies:
RUN pip install cython==0.29.19 
RUN pip install pysam==0.15.4 biopython==1.76 numpy==1.16.6 mailer==0.8.1

RUN conda config --add channels r
RUN conda config --add channels bioconda
#Needed for HiddenMarkov:
RUN conda config --add channels conda-forge

RUN conda upgrade conda

#Use R 3.5.1 because latest R is 4.0.0 and doesn't have updated packages
#Use bwa 0.7.17-0:
RUN conda install bwa=0.7.17-0 r=3.5.1 r-r.oo r-r.methodss3 r-hiddenmarkov

COPY ./ /msg
COPY ./cmdline /msg/cmdline
COPY ./dependencies/samtools-0.1.9.tar.bz2 /usr/opt
#An older Stampy (1.0.23) is packaged with older commits as Stampy-latest.tgz
# but I'll test with the newest Stampy (1.0.32) for now:
COPY ./dependencies/stampy-1.0.32.tgz /usr/opt

WORKDIR /usr/opt

#Might not need this if we can get the MSG makefile working with Docker:
# install samtools
RUN tar -jvxf samtools-0.1.9.tar.bz2 && \ 
    cd samtools-0.1.9 && \
    make

#Might not need this if we can get the MSG makefile working with Docker:
ENV PATH /usr/opt/samtools-0.1.9:$PATH

#Might not need this if we can get the MSG makefile working with Docker, and
# add Stampy installation to the makefile:
RUN tar -zvxf stampy-1.0.32.tgz && \
    cd stampy-1.0.32 && \
    make

ENV PATH /usr/opt/stampy-1.0.32:$PATH

WORKDIR /msg

RUN gcc -O2 -lm -o countalleles countalleles.c && \
    gcc -O2 -lm -o hmmprobs hmmprobs.c 

RUN sh test_dependencies_R.sh

#It would be nice to be able to use the existing makefile instead of
# recapitulating parts of it with RUN statements:
#RUN make
#ENV PATH /msg/samtools-0.1.9:$PATH

WORKDIR /data
