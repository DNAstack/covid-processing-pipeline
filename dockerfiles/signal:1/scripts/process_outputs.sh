#!/bin/bash

#!/bin/bash

trap "rm new_samplename.txt header.txt" EXIT

usage() {
cat << EOF

Usage: $0 -s samplename -m min_frequency -v reference_basename

  -h    Print this message and exit
  -s    Samplename
  -m    High confidence minimum frequency threshold
  -v 	Viral reference basename

EOF
}

while getopts "hs:m:v:" OPTION; do
	case $OPTION in
		h) usage; exit;;
		s) SAMPLENAME=$OPTARG;;
		m) MIN_FREQ_HIGH_CONFIDENCE_THRESHOLD=$OPTARG;;
		v) VIRAL_REFERENCE_BASE=$OPTARG;;
		\?) usage; exit;;
	esac
done

if [ -z "${SAMPLENAME}" ] || [ -z "${MIN_FREQ_HIGH_CONFIDENCE_THRESHOLD}" ] || [ -z "${VIRAL_REFERENCE_BASE}" ]; then
	usage
	exit 1
fi


mv -f output/* .

cp "${SAMPLENAME}/core/${SAMPLENAME}_ivar_variants.tsv" .
cp "${SAMPLENAME}/core/${SAMPLENAME}.consensus.fa" .
cp "${SAMPLENAME}/core/${SAMPLENAME}_viral_reference.bam" "${SAMPLENAME}.bam"
samtools index "${SAMPLENAME}.bam"

ivar_variants_to_vcf.py \
	"${SAMPLENAME}_ivar_variants.tsv" \
	"${SAMPLENAME}.vcf"

chrom_line=$(bcftools view -h "${SAMPLENAME}.vcf" | tail -1)
echo "$(echo "${chrom_line}" | awk '{print $10}') ${SAMPLENAME}" > new_samplename.txt

bcftools view --no-version -h "${SAMPLENAME}.vcf" | sed '$ d' > header.txt
echo -e "##processing_pipeline=https://github.com/jaleezyy/covid-19-signal/tree/$SIGNAL_VERSION\n$chrom_line" >> header.txt
bcftools reheader \
	-h header.txt \
	-s new_samplename.txt \
	"${SAMPLENAME}.vcf" \
| bgzip > "${SAMPLENAME}.low_freq.vcf.gz"
tabix "${SAMPLENAME}.low_freq.vcf.gz"

bcftools view \
	--no-version \
	-i "FORMAT/ALT_FREQ > ${MIN_FREQ_HIGH_CONFIDENCE_THRESHOLD}" \
	-f "PASS" \
	"${SAMPLENAME}.low_freq.vcf.gz" \
| bgzip > "${SAMPLENAME}.vcf.gz"
tabix "${SAMPLENAME}.vcf.gz"

new_fasta_header=">${SAMPLENAME} $(head -1 "${SAMPLENAME}.consensus.fa" | tr -d ">") ${VIRAL_REFERENCE_BASE}"
sed -i "1s/^.*$/$new_fasta_header/" "${SAMPLENAME}.consensus.fa"
