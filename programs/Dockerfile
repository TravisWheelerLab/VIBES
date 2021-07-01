FROM ubuntu:20.04

# install necessary packages, tools. noninteractive stops tzdata from asking for
# geographic time zone (stops Docker)
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install \
        build-essential \
        curl \
        perl \
        python3 \
        python3-matplotlib \
        wget
RUN cpan install CPAN && \
    cpan reload cpan && \
    cpan install Moose

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

RUN mkdir -p /programs
ADD . /programs/

WORKDIR /

