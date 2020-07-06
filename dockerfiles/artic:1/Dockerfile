# TAG dnastack/artic:1
FROM continuumio/miniconda3:4.7.12

MAINTAINER Heather Ward <heather@dnastack.com>

# Workaround for java install errors
RUN mkdir -p /usr/share/man/man1

RUN apt-get -qq update && \
	apt-get -qq install \
	openjdk-11-jdk \
	curl \
	bcftools \
	tabix \
	samtools

# Install nextflow
RUN curl -fsSL get.nextflow.io | bash && \
	mv nextflow /bin/

# Get connor-lab workflow
RUN git clone https://github.com/connor-lab/ncov2019-artic-nf.git /opt/ncov2019-artic-nf

WORKDIR /opt/ncov2019-artic-nf

ENV ARTIC_VERSION 9a91994b25adcb6b1ac6b5af48321fad37dcce2a
RUN git checkout ${ARTIC_VERSION}

# Generate conda cache using dummy file
COPY ./ERR4007730_1.400.fastq /mnt/
RUN cd /mnt && \
	nextflow run /opt/ncov2019-artic-nf/main.nf -profile conda --cache /opt/conda/envs/ --medaka --prefix ERR4007730 --basecalled_fastq /mnt && \
	rm -r /mnt/*

RUN git clone https://github.com/artic-network/artic-ncov2019.git /opt/artic-ncov2019 && \
	samtools faidx /opt/artic-ncov2019/primer_schemes/nCoV-2019/V2/nCoV-2019.reference.fasta

# Replace the artic.nf file with our version; prevents downloading of the scheme repo
COPY ./artic.nf /opt/ncov2019-artic-nf/modules/artic.nf

COPY ./scripts/* /usr/local/bin/

