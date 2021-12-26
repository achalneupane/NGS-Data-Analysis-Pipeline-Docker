#### Building DOCKER image for GATK-4.0 / baserecalibrator
# Specify the base image -- here we're using one that bundles the OpenJDK version of Java 8 on top of a generic Debian Linux OS
# OpenJDK is an open-source implementation of the Java Platform, Standard Edition

# Specify the base image 
FROM openjdk:8-jre
##Alternatives: 
## FROM ubuntu:xenial  # as in the mgibio/alignment
## FROM phusion/baseimage:0.11 as builder 

LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.schema-version="1.01"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/baserecalibrator"
LABEL org.label-schema.description="Image for running GATKv4/baserecalibrator at the MGI"

#Set the working directory to be used when the docker gets run
WORKDIR /tmp   

# Do a few updates of the base system 
RUN apt-get update -y && \
	apt-get install -y \
    ant \
    apt-utils \
    build-essential \
    bzip2 \
    default-jdk \
    default-jre \
    gcc-multilib \
    git \
    libncurses5-dev \
    libnss-sss \
    libz-dev \
	libbz2-dev \
	liblzma-dev \
	nodejs \
    python-dev \
    python-pip \
    r-base \
    time \
	tzdata \
    unzip \
    wget \
    zlib1g-dev
	
RUN apt-get update -y \
	&& apt-get clean all \
	&& rm  -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/*
	
# Define a timezone so Java works properly
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime \
    && echo "America/Chicago" > /etc/timezone \
    && dpkg-reconfigure --frontend noninteractive tzdata

# Install R Packages
RUN Rscript -e 'install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)'

# Install the ggplot2 library and a few other dependencies we want to have available
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"


##################
#### GATK 4.1.2.0 ####
##################

#Install GATK4
RUN cd /tmp/ \
	&& wget -q 'https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip' \
    && unzip gatk-4.1.2.0.zip \
	&& cp -av /tmp/gatk-4.1.2.0 /usr/bin/

## Export gatk directory ti environmental variables
## if using this export then later call GATK using  ${GATKv4}
ENV GATKv4 /usr/bin/gatk-4.1.2.0/gatk  
## if using this export then just call GATK using ./gatk from your WORKDIR
#ENV PATH=/usr/bin/gatk-4.1.2.0:$PATH   

# helper scripts
COPY ./entrypoint.sh /usr/bin/entrypoint.sh
ENTRYPOINT ["/usr/bin/entrypoint.sh"]
