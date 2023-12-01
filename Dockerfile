FROM ubuntu:20.04

# install necessary packages, tools. noninteractive stops tzdata from asking for
# geographic time zone (stops Docker)
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install \
        build-essential \
        curl \
        perl \
        python3-pip \
        wget \
        git \
        autoconf

# install pandas
RUN python3 -m pip install pandas numpy

ENV PATH "$PATH:/tool_bins"

# install HMMER
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar zxf hmmer.tar.gz && \
    cd hmmer-3.4 && \
    ./configure --prefix /tool_bins && \
    make && \
    make check && \
    make install && \
    cd easel && make install

ENV PATH "$PATH:/tool_bins/bin"

# install BATH
RUN git clone https://github.com/TravisWheelerLab/BATH && \
   cd BATH && \
   git clone https://github.com/TravisWheelerLab/easel && \
   cd easel && \
   git checkout BATH && \
   cd ../ && \
   autoconf && \
   ./configure && \
   make && \
   # make check && \
   make install

ENV PATH "$PATH:/BATH/src/"

# Install AWS command line tools
RUN apt-get update && apt-get -y install unzip
WORKDIR /opt
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN aws --version

# copy workflow scripts and add to $PATH
RUN mkdir vibes_scripts/
COPY ../nextflow_workflow/bin/ vibes_scripts/
ENV PATH "$PATH:/vibes_scripts"

WORKDIR /
