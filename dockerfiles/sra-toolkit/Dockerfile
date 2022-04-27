# TAG dnastack/sra-toolkit:2.10.7
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		wget

ENV SRATOOLKIT_VERSION 2.10.7
RUN wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz" && \
	tar -zxvf "sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz" --directory /opt/ && \
	rm "sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz"
ENV PATH $PATH:/opt/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64/bin/

# Generate vdb-config; dummy GUID
RUN vdb-config --restore-defaults && \
	echo '/LIBS/GUID = "a924eb97-2c35-47ed-89da-5b82f13e6804"' >> /root/.ncbi/user-settings.mkfg
