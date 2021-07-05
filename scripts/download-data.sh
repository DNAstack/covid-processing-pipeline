#!/bin/bash


function usage() {
cat << EOF

Download files required to run the COVID-19 analysis pipeline
Usage: $0 [options]

  -h    Print this message and exit
  -r    Redownload files (replace existing files if they exist). Default is to skip downloading existing files
  -f		Download data needed for the 'from_fastq' workflow
  -a		Download data needed for the 'from_assmebly' workflow

EOF
}

function getFile() {
FILE_URL=$1
OUTPUT_FILE_NAME="${DATA_DIR}/$(basename "${FILE_URL}")"

if [ -e "${OUTPUT_FILE_NAME}" ] && [ "${REDOWNLOAD}" = "TRUE" ] || [ ! -e "${OUTPUT_FILE_NAME}" ]; then
	retries=3
	while [ "${retries}" -gt 0 ]; do
		curl -X GET "${FILE_URL}" --output "${OUTPUT_FILE_NAME}" && break
		((retries--))
	done
fi
}


## Config
DATA_DIR=data

### From fastq data
KRAKEN2_DB=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/Kraken2.tar.gz
COMPOSITE_REFERENCE=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/composite_human_viral_reference.tar.gz
VIRAL_REFERENCE=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.fasta
VIRAL_REFERENCE_FEATURE_COORDS=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.gff3

### From assembly data
VIRAL_REFERENCE=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.fasta
SNPEFF_DB="https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/NC_045512.2.zip"


while getopts "hrfa" OPTION; do
	case $OPTION in
		h) usage; exit;;
		r) REDOWNLOAD=TRUE;;
		f) DOWNLOAD_FROM_FASTQ_DATA=TRUE;;
		a) DOWNLOAD_FROM_ASSEMBLY_DATA=TRUE;;
		\?) usage; exit;;
	esac
done


## Main
if [ "${DOWNLOAD_FROM_FASTQ_DATA}" == "TRUE" ]; then
	getFile "${KRAKEN2_DB}"
	getFile "${COMPOSITE_REFERENCE}"
	getFile "${VIRAL_REFERENCE}"
	getFile "${VIRAL_REFERENCE_FEATURE_COORDS}"
fi

if [ "${DOWNLOAD_FROM_ASSEMBLY_DATA}" == "TRUE" ]; then
	getFile "${VIRAL_REFERENCE}"
	getFile "${SNPEFF_DB}"
fi
