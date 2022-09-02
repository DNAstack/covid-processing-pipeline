#!/bin/bash

set -eEo pipefail

IFS=$'\n'

usage() {
cat << EOF

Usage: $0 -a accession -p artic_primer_version

OPTIONS
  -h    Display this message and exit
  -a    Sample accession
  -p    Inferred ARTIC primer version

EOF
}

while getopts "ha:p:" OPTION; do
	case $OPTION in
		h) usage; exit 1;;
		a) ACCESSION=$OPTARG;;
		p) ARTIC_PRIMER_VERSION=$OPTARG;;
		\?) usage; exit 1;;
	esac
done

if [ -z "${ACCESSION}" ]; then
	usage
	echo "[ERROR] Must provide sample accession"
	exit 1
fi

METADATA=/tmp/${ACCESSION}.runinfo.csv
META_PARSED=/tmp/meta.parsed.csv

## Get run metadata
count=0
retc=1
echo "Attempting metadata access"

set +e
while [[ "$count" -lt 3  && ( "$retc" -ne 0 || -z "${runinfo}" ) ]]; do
	sleep 60
	echo "[$(date)] Metadata: Attempt [$count]"
	runinfo=$(esearch -db sra -query "${ACCESSION}" | efetch -format runinfo | sed '/^$/,+1 d')
	retc=$?

	echo "retc is [$retc]"
	echo
	echo "runinfo is [$runinfo]"
	echo

	((++count))
done
set -e
if [ "$retc" -ne 0 ] || [ -z "$runinfo" ]; then
	>&2 echo "[ERROR] Could not retrieve runinfo"
	exit 1
else
	echo "${runinfo}" > "${METADATA}"
fi

echo "Metadata successfully downloaded"
cat "${METADATA}" | sed -r ':r; s/("[^",]+),([^",]*)/\1 \2/g; tr; s/"//g' > "${META_PARSED}"

## Get BioSample metadata
biosample_column=$(head -1 "${META_PARSED}" | tr ',' '\n' | grep -n BioSample | cut -d : -f 1)
if [ -n "${biosample_column}" ]; then
	biosample=$(tail -1 "${META_PARSED}" | cut -d , -f ${biosample_column})
	if [ -n "${biosample}" ]; then
		count=0
		retc=1
		echo "Attempting attribute access"

		set +e
		while [ "$count" -lt 3 ] && [ "$retc" -ne 0 ]; do
			sleep 5
			echo "Attribute access: attempt [$count]"
			attributes=$(esearch -db biosample -q "${biosample}" | efetch -mode xml -json | jq -r '.BioSampleSet.BioSample.Attributes.Attribute[] | "\(.harmonized_name)\t\(.attribute_name)\t\(.content)"')
			retc=$?

			((++count))
		done
		set -e

		echo "Attributes successfully downloaded"

		if [ $retc -ne 0 ]; then
			echo "[WARN] Failed to retrieve BioSample data"
		else
			set +e
			# use only official BioSample attributes
			official_attributes=$(echo -e "${attributes}" | awk '{if (! ($1 == "null")) print}' | grep -i -v -f "${MISSING_DATA_PATTERNS}")
			set -e

			if [ "$(echo "${official_attributes}" | wc -c)" -gt 1 ]; then
				# prefer host scientific name
				if [ "$(echo "${official_attributes}" | grep -w -c "^host")" -gt 1 ]; then
					official_attributes=$(echo "${official_attributes}" | grep -v "host scientific name")
				fi

				# confirm no duplicate attributes
				if [ "$(echo "${official_attributes}" | cut -f 1 | sort | uniq -c | awk '{print $1}' | sort | uniq | wc -l)" -gt 1 ]; then
					>&2 echo "[ERROR] duplicate attribute in extended metadata found"
					>&2 echo "Official attributes:"
					>&2 echo "--------------------"
					>&2 echo "${official_attributes}"
					>&2 echo "--------------------"
					exit 1
				fi

				for attribute in $(echo -e "${official_attributes}"); do
					att_name=$(echo "${attribute}" | cut -f 1)
					# quote fields with commas
					att_val=$(echo "${attribute}" | cut -f 3 | sed 's/\(.*,.*\)/"\1"/')
					header="${header},${att_name}"
					meta_line="${meta_line},${att_val}"
				done

				echo -e "${header}\n${meta_line}" | sed -e "s#^,##g" -e '1s# #_#g' -e '1s#\(.*\)#\L\1#g' \
					> additional_meta.csv

				if [[ -n "${ARTIC_PRIMER_VERSION}" ]]; then
					paste -d , "${METADATA}" additional_meta.csv <(echo -e "inferred_artic_primer_version\n${ARTIC_PRIMER_VERSION}") > "${ACCESSION}.meta.csv"
				else
					paste -d , "${METADATA}" additional_meta.csv > "${ACCESSION}.meta.csv"
				fi
				exit 0
			else
				echo "[WARN] No official attributes found"
			fi
		fi
	else
		echo "[WARN] BioSample column empty"
	fi
else
	echo "[WARN] No BioSample column detected"
fi

# If we failed to access additional attributes using BioSample
if [[ -n "${ARTIC_PRIMER_VERSION}" ]]; then
	paste -d , "${METADATA}" <(echo -e "inferred_artic_primer_version\n${ARTIC_PRIMER_VERSION}") > "${ACCESSION}.meta.csv"
else
	mv "${METADATA}" "${ACCESSION}.meta.csv"
fi
