FROM openjdk:8-jre

LABEL maintainer="Achal Neupane <achalneupane@wustl.edu>"
LABEL org.label-schema.schema-version="1.0"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/bbduk"
LABEL org.label-schema.description="Image for creating per lane fastq"


RUN    apt-get update \
    && apt-get install -y \
        build-essential \
        openjdk-11-jre \
        wget \
    && rm -rf /var/lib/apt/lists/*

ENV BBMAP_VERSION 38.87
ENV BBMAP_DIR /usr/local/bbmap

WORKDIR /usr/local
RUN set -eux; \
       wget -O bbmap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download \
    && tar -zxf bbmap.tar.gz \
    && rm bbmap.tar.gz

ENV PATH ${BBMAP_DIR}:${PATH}

## Add ENV for Shell scripts from BBMAP
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
# CMD ["/bin/bash"]