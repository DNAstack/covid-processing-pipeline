#!/bin/bash


function usage() {
cat << EOF

Download files required to run the COVID-19 analysis pipeline
Usage: $0 [options]

  -h    Print this message and exit
  -r    Redownload files (replace existing files if they exist). Default is to skip downloading existing files

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
KRAKEN2_DB=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/Kraken2.tar.gz
HUMAN_REFERENCE=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.tar.gz
VIRAL_REFERENCE=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.fasta
VIRAL_REFERENCE_FEATURE_COORDS=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.gff3
# For breseq
VIRAL_GBK=https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.gbk


while getopts "hr" OPTION; do
	case $OPTION in
		h) usage; exit;;
		r) REDOWNLOAD=TRUE;;
		\?) usage; exit;;		
	esac
done


## Main

getFile "${KRAKEN2_DB}" 
getFile "${HUMAN_REFERENCE}"
getFile "${VIRAL_REFERENCE}"
getFile "${VIRAL_REFERENCE_FEATURE_COORDS}"
getFile "${VIRAL_GBK}"

