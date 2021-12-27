FROM openjdk:8-jre
LABEL maintainer="Achal Neupane <achalneupane@wustl.edu>"

LABEL org.label-schema.schema-version="1.0"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/bbduk"
LABEL org.label-schema.description="Image for creating interleaved fastq"

RUN apt-get update -y \
&& apt-get install --no-install-recommends -y \
 make \
 gcc \
 g++ \
 libz-dev \
 libbz2-dev \
 liblzma-dev \
 ncurses-dev \
 bc \
 libnss-sss \
 time 

ENV BBMAP_URL http://downloads.sourceforge.net/project/bbmap/BBMap_38.87.tar.gz
RUN cd /usr/local/ \
    && wget -O bbmap.tar.gz ${BBMAP_URL} \
    && tar -xvzf bbmap.tar.gz \
    && rm bbmap.tar.gz
ENV PATH="/usr/local/bbmap:${PATH}"

CMD ["bbduk.sh"]