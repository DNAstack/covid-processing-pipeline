# TAG dnastack/signal:1
FROM continuumio/miniconda3:4.7.12

# Based on https://github.com/jaleezyy/covid-19-signal/blob/master/resources/Dockerfile_pipeline
MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		curl \
		wget \
		git \
		build-essential \
		bedtools \
		bcftools \
		tabix \
		samtools

RUN conda create --name snakemake --channel conda-forge --channel bioconda snakemake=5.11.2 pandas

RUN git clone https://github.com/jaleezyy/covid-19-signal.git /opt/covid-19-signal

WORKDIR /opt/covid-19-signal

ENV SIGNAL_VERSION f7c6f65e32fb1c60873b361e6ca0c8c7241d1a44
RUN git checkout ${SIGNAL_VERSION}

# Force creation of cached conda environments
RUN cp example_config.yaml config.yaml && \
	cp example_sample_table.csv sample_table.csv && \
	mkdir data && \
	touch data/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz \
		data/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz \
		data/MN908947.3.fasta \
		data/MN908947.3.gff3 && \
	conda run \
		-n snakemake \
		snakemake \
		--verbose \
		--use-conda \
		--create-envs-only \
		--cores 12 \
		--conda-prefix $HOME/.snakemake \
		-s Snakefile all && \
	rm -rf data config.yaml sample_table.csv

COPY ./scripts/* /usr/local/bin/
