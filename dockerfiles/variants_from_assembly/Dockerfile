# TAG dnastack/variants_from_assembly:latest
FROM ubuntu:bionic

RUN apt-get -qq update && \
	apt-get -qq install \
	curl \
	wget \
	unzip \
	ncbi-blast+ \
	snp-sites \
	bcftools \
	vcftools \
	openjdk-8-jre \
	gawk

# mview 1.67
RUN wget https://github.com/desmid/mview/archive/v1.67.tar.gz && \
    tar -zxvf v1.67.tar.gz --directory /opt/ && \
    rm v1.67.tar.gz
ENV PERL5LIB=$PERL5LIB:/opt/mview-1.67/lib/
RUN ln -s /opt/mview-1.67/bin/mview /usr/local/bin/

# snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
	unzip -d /opt/ snpEff_latest_core.zip && \
	rm snpEff_latest_core.zip

COPY ./scripts/* /usr/local/bin/
