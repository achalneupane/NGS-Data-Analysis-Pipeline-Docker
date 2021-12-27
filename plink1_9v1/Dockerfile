FROM openjdk:8-jre

LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.name="achalneupane/plink1.9v1"
LABEL org.label-schema.description="Image for plink"


RUN apt-get update
COPY plink1.9 /usr/bin/

ENV PATH=/usr/bin/:$PATH

## To run
# docker run -ti achalneupane/plink1.9v1 plink1.9 --help
