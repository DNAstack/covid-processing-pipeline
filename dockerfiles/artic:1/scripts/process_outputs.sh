#!/bin/bash

trap "rm new_samplename.txt header.txt" EXIT

usage() {
cat << EOF

Usage: $0 -s samplename

  -h    Print this message and exit
  -s    Samplename

EOF
}

while getopts "hs:" OPTION; do
	case $OPTION in
		h) usage; exit;;
		s) SAMPLENAME=$OPTARG;;
		\?) usage; exit;;
	esac
done

if [ -z "${SAMPLENAME}" ]; then
	usage
	echo "[ERROR]: Provide samplename"
	exit 1
fi

# High confidence variants
# Change samplename from SAMPLE to the run ID
chrom_line=$(bcftools view -h ${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.pass.vcf.gz | tail -1)
echo "SAMPLE ${SAMPLENAME}" > new_samplename.txt

bcftools view --no-version -h ${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.pass.vcf.gz | sed '$ d' > header.txt
echo -e "##processing_pipeline=https://github.com/connor-lab/ncov2019-artic-nf/tree/$ARTIC_VERSION\n$chrom_line" >> header.txt
bcftools reheader \
	-h header.txt \
	-s new_samplename.txt \
	${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.pass.vcf.gz \
	-o ${SAMPLENAME}.vcf.gz
tabix ${SAMPLENAME}.vcf.gz


# Medaka variants
# Change samplename from SAMPLE to the run ID
chrom_line=$(bcftools view -h ${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.merged.vcf.gz | tail -1)

bcftools view --no-version -h ${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.merged.vcf.gz | sed '$ d' > header.txt
echo -e "##processing_pipeline=https://github.com/connor-lab/ncov2019-artic-nf/tree/$ARTIC_VERSION\n$chrom_line" >> header.txt
bcftools reheader \
	-h header.txt \
	-s new_samplename.txt \
	${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.merged.vcf.gz \
	-o ${SAMPLENAME}.medaka.vcf.gz
tabix ${SAMPLENAME}.medaka.vcf.gz


cp ${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}*.consensus.fasta ${SAMPLENAME}.consensus.fa
sed -i '1!b;s/_execution//' ${SAMPLENAME}.consensus.fa

find -regex ".*/${SAMPLENAME}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/${SAMPLENAME}[^.]*.sorted.bam" -regextype sed -exec cp {} ./${SAMPLENAME}.bam \;
samtools index ${SAMPLENAME}.bam
