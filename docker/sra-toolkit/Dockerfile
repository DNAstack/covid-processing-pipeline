# TAG dnastack/sra-toolkit:2.10.7
FROM ubuntu:bionic

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		tabix \
		jq \
		ncbi-entrez-direct \
		rs

ENV SRATOOLKIT_VERSION 2.10.7
RUN wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz" && \
	tar -zxvf "sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz" --directory /opt/ && \
	rm "sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz"
ENV PATH $PATH:/opt/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64/bin/

RUN wget http://ftp.debian.org/debian/pool/main/n/ncbi-entrez-direct/ncbi-entrez-direct_10.9.20190219+ds-1+b10_amd64.deb && \
	apt-get -qq install ./ncbi-entrez-direct_10.9.20190219+ds-1+b10_amd64.deb && \
	rm ./ncbi-entrez-direct_10.9.20190219+ds-1+b10_amd64.deb

RUN mkdir -p /root/.ncbi && \
	echo '/LIBS/GUID = "a924eb97-2c35-47ed-89da-5b82f13e6804"' >> /root/.ncbi/user-settings.mkfg

COPY ./scripts/* /usr/local/bin/

ENV MISSING_DATA_PATTERNS /opt/missing_data_patterns.txt
COPY ./missing_data_patterns.txt "${MISSING_DATA_PATTERNS}"
