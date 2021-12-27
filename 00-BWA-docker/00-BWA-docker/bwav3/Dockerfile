FROM openjdk:8-jre
LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.schema-version="2.2"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/bwaref2"
LABEL org.label-schema.description="Image for bwa"

ENV BWA_VERSION 0.7.17
ENV PICARD_VERSION 2.18.21
ENV SAMTOOLS_VERSION 1.9

WORKDIR /tmp

RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       make \
       gcc \
       g++ \
       libz-dev \
       libbz2-dev \
       liblzma-dev \
       ncurses-dev \
       libnss-sss \
       time \
  && wget -q https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
  && tar xjvf bwa-${BWA_VERSION}.tar.bz2 \
  && cd /tmp/bwa-${BWA_VERSION}/ \
  && make \
  && cp -av /tmp/bwa-${BWA_VERSION}/bwa /usr/bin/ \
  && cd /tmp \
  && wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd /tmp/samtools-${SAMTOOLS_VERSION}/ \
  && make \
  && cp -av /tmp/samtools-${SAMTOOLS_VERSION}/samtools /usr/bin/ \
  && wget -q -O /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
  && ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime \
  && echo "America/Chicago" > /etc/timezone \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/* 

# This makes the image crazy large -- will find a workaround
# COPY human_g1k_v37_decoy* /usr/local/

COPY ./entrypoint.sh /usr/local/bin/
ENV PICARD /usr/bin/picard.jar

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
# CMD ["/bin/bash"]