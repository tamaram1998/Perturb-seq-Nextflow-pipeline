# Use a base image with Ubuntu 22.04
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-11-jdk \
    wget \
    unzip \
    python3-pip \
    build-essential \
    curl \
    vim-common \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Install FastQC
RUN apt-get update && \
    apt-get install -y fastqc

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.zip && \
    unzip 2.7.11b.zip && \
    cd STAR-2.7.11b/source && \
    sed -i 's/-mavx2//g' Makefile && \
    make STAR && \
    mv STAR /usr/local/bin/ && \
    cd ../../ && \
    rm -rf STAR-2.7.11b*

# Install featureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.7/subread-2.0.7-Linux-x86_64.tar.gz && \
    tar -xzf subread-2.0.7-Linux-x86_64.tar.gz && \
    mv subread-2.0.7-Linux-x86_64/bin/* /usr/local/bin/ && \
    rm -rf subread-2.0.7-Linux-x86_64*

# Install MultiQC
RUN pip3 install multiqc

# Set the working directory
WORKDIR /data

# Copy files
COPY . .

# Specify default command
CMD ["nextflow", "run", "/pipeline/main.nf"]
