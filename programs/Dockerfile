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

# install commander
RUN mkdir tool_bins/ && cd tool_bins/ && curl -L https://github.com/commander-cli/commander/releases/download/v2.3.0/commander-linux-amd64 -o commander && chmod +x commander
ENV PATH "$PATH:/tool_bins"

# install HMMER
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar zxf hmmer.tar.gz && \
    cd hmmer-3.3.2 && \
    ./configure --prefix /tool_bins && \
    make && \
    make check && \
    make install && \
    cd easel && make install

ENV PATH "$PATH:/tool_bins/bin"

# install frahmmer
RUN git clone https://github.com/TravisWheelerLab/FraHMMER.git && \
   cd FraHMMER && \
   git clone https://github.com/TravisWheelerLab/easel && \
   (cd easel; git checkout develop) && \
   autoconf && \
   ./configure --prefix /tool_bins && \
   make && \
   make check && \
   make install

# Install AWS command line tools
RUN apt-get update && apt-get -y install unzip
WORKDIR /opt
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN aws --version

RUN mkdir -p /programs
COPY python /programs/python
COPY perl /programs/perl

WORKDIR /
