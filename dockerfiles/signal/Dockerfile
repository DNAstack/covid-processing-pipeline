# TAG dnastack/signal:09372f8
FROM continuumio/miniconda3:4.9.2

# Based on https://github.com/jaleezyy/covid-19-signal/blob/master/resources/Dockerfile_pipeline
MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		curl \
		wget \
		git \
		build-essential \
		libtbb2 \
		tabix \
		bcftools \
		zip && \
	conda init bash

# install snakemake and system dependencies
RUN conda create --name snakemake --channel conda-forge --channel bioconda snakemake==6.0.5 pandas

# clone the workflow repo (tag later to specific release to match software version)
RUN git clone https://github.com/jaleezyy/covid-19-signal

# set workdir to pipeline folder in repo
WORKDIR /covid-19-signal

ENV SIGNAL_VERSION 09372f8e699adacc7c3e41844606d117bd65b952
RUN git checkout ${SIGNAL_VERSION}

# install all snakemake rule conda envs using dummy snakemake file
COPY ./dependencies resources/dependencies
RUN conda run \
	-n snakemake snakemake \
	--verbose \
	--use-conda \
	--conda-prefix=$HOME/.snakemake \
	-s resources/dependencies \
	--cores 8

COPY ./scripts/* /usr/local/bin/

# Do not automatically update pangolin; can both break output and reach github rate limits
RUN sed -i 's/\(^[^d]*update_pangolin()\)/#\1/' /covid-19-signal/scripts/assign_lineages.py
