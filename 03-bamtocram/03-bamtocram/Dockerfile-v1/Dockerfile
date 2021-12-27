FROM phusion/baseimage:0.11 as builder
ENV SAMTOOLS_VERSION 1.9

#set the working directory to be used when the docker gets run
WORKDIR /tmp

# Do a few updates of the base system
RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       wget \
       make \
       g++ \
       ncurses-dev \
       libz-dev \
       libbz2-dev \
       liblzma-dev \
  && wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd /tmp/samtools-${SAMTOOLS_VERSION}/ \
  && make
  
  
FROM phusion/baseimage:0.11
LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.schema-version="1.1"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/bamtocram"
LABEL org.label-schema.description="Image for samtools"

RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
 	   libnss-sss \
       time \
       tzdata \ 
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/* 

ENV SAMTOOLS_VERSION 1.9
ENV TZ "America/Chicago"

COPY --from=builder /tmp/samtools-${SAMTOOLS_VERSION}/samtools  /usr/bin/
COPY ./entrypoint.sh /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]