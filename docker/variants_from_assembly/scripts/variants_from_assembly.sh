#!/bin/bash

function usage() {
cat << EOF

Usage: $0 -s sample_id -a sample_assembly -r reference_genome -i reference_genome_id

EOF
}

function splitFasta() {
	MULTIFASTA=$1
	OUTPUT_DIR=$2

	awk -v OUTPUT_DIR="${OUTPUT_DIR}" '{if (substr($0, 1, 1)==">") {num_digits=index($0, " ")-2;filename=gensub(/\//, "_", "g", (substr($0,2,num_digits) ".fa"))} print $0 > OUTPUT_DIR"/"filename}' \
		< "${MULTIFASTA}"
}

while getopts "hs:a:r:i:" OPTION; do
	case $OPTION in
		h) usage; exit;;
		s) SAMPLE=$OPTARG;;
		a) ASSEMBLY=$OPTARG;;
		r) REFERENCE=$OPTARG;;
		i) REFERENCE_GENOME_ID=$OPTARG;;
		\?) usage; exit;;
	esac
done

if [ -z "${SAMPLE}" ] || [ -z "${ASSEMBLY}" ] || [ -z "${REFERENCE}" ] || [ -z "${REFERENCE_GENOME_ID}" ]; then
	usage; exit 1
fi

# Change assembly header if it does not match samplename
# This will only impact the ability to create the variants file, not the output assembly
assembly_samplename=$(head -1 "${ASSEMBLY}" | cut -d ' ' -f 1 | sed 's/^>//')
if [ ! "${assembly_samplename}" = "${SAMPLE}" ]; then
	sed -i "s~^>${assembly_samplename}~>${SAMPLE}~" "${ASSEMBLY}"
fi

## VCF header lines
VCF_GENERIC_HEADER="##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##contig=<ID=NC_045512,length=29903>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
VCF_REFERENCE_LINE="##reference=file:///${REFERENCE_GENOME_ID}.fasta"
VCF_BLASTN_VERSION="##blastnVersion=\"2.6.0+\""
VCF_MVIEW_VERSION="##mviewVersion=\"1.67\""
VCF_SNP_SITES_VERSION="##snp-sitesVersion=\"2.3.2\""
VCF_BLASTN_COMMANDLINE="##blastnCmd=\"blastn -query ${REFERENCE_GENOME_ID}.fasta -subject ${SAMPLE}.fasta -outfmt 0 | mview -in blast -out fasta > ${SAMPLE}.aln.fasta\""
VCF_SNP_SITES_COMMANDLINE="##snp-sitesCmd=\"snp-sites ${SAMPLE}.aln.fasta -v -o ${SAMPLE}.vcf\""


mkdir tmp

# Allow renaming the vcf chr from 1 to REFERENCE_GENOME_ID
echo -e "1\t${REFERENCE_GENOME_ID}" > tmp/chr_map_file.txt

# Including the full reference sequence in the BLAST subject file forces the output positions to match up to the base numbering of the full
# reference - otherwise numbering is based on the first position that the subject matches in the reference, which could be > position 1 of
# the reference (position in VCF would be wrong)
cat ${ASSEMBLY} ${REFERENCE} > ${SAMPLE}.blast_subject.fasta

blastn \
	-query ${REFERENCE} \
	-subject ${SAMPLE}.blast_subject.fasta \
	-outfmt 0 | \
mview \
	-in blast \
	-out fasta | \
sed -e "s/[0-9]://g" \
	-e "s/^>${SAMPLE}[^ ]*/>${SAMPLE}/" > \
	${SAMPLE}.aligned.fasta

# Get rid of the duplicate reference sequences in the alignment
splitFasta "${SAMPLE}.aligned.fasta" ./tmp

if [ -e "tmp/${SAMPLE}.fa" ]; then
	# a passing alignment exists
	# snp-sites requires that the reference sequence is first
	cat ${REFERENCE} "tmp/${SAMPLE}.fa" > "tmp/${SAMPLE}.aln.fasta"

	snp-sites \
		"tmp/${SAMPLE}.aln.fasta" \
		-v \
		-o "tmp/${SAMPLE}.vcf"

	if [ -e "tmp/${SAMPLE}.vcf" ]; then
		# at least one variant exists for this sample
		bcftools view \
			--no-version \
			--no-update \
			-s "${SAMPLE}" \
			"tmp/${SAMPLE}.vcf" | \
		bcftools annotate \
			--no-version \
			--rename-chrs tmp/chr_map_file.txt \
			- \
			-o "tmp/${SAMPLE}.renamed.vcf"

		# Clean up header; add tool versions
		header=$(bcftools view \
			--no-version \
			--no-update \
			-h \
			"tmp/${SAMPLE}.renamed.vcf" | head -n -1)
		CHROM_LINE=$(bcftools view \
			--no-version \
			--no-update \
			-h \
			"tmp/${SAMPLE}.renamed.vcf" | tail -1)

		echo -e "${header}\n${VCF_REFERENCE_LINE}\n${VCF_BLASTN_VERSION}\n${VCF_MVIEW_VERSION}\n${VCF_SNP_SITES_VERSION}\n${VCF_BLASTN_COMMANDLINE}\n${VCF_SNP_SITES_COMMANDLINE}\n${CHROM_LINE}" > "tmp/${SAMPLE}.header.txt"
		bcftools reheader \
			-h "tmp/${SAMPLE}.header.txt" \
			"tmp/${SAMPLE}.renamed.vcf" \
			-o "${SAMPLE}.vcf"

	else
		# no variants exist for this sample
		CHROM_LINE="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SAMPLE}"
		echo -e "${VCF_GENERIC_HEADER}\n${VCF_REFERENCE_LINE}\n${VCF_BLASTN_VERSION}\n${VCF_MVIEW_VERSION}\n${VCF_SNP_SITES_VERSION}\n${VCF_BLASTN_COMMANDLINE}\n${VCF_SNP_SITES_COMMANDLINE}\n${CHROM_LINE}" \
				> "${SAMPLE}.vcf"

	fi
else
	# no passing alignments for this sample
	CHROM_LINE="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SAMPLE}"
	echo -e "${VCF_GENERIC_HEADER}\n${VCF_REFERENCE_LINE}\n${VCF_BLASTN_VERSION}\n${VCF_MVIEW_VERSION}\n${VCF_SNP_SITES_VERSION}\n${VCF_BLASTN_COMMANDLINE}\n${VCF_SNP_SITES_COMMANDLINE}\n${CHROM_LINE}" \
			> "${SAMPLE}.vcf"

fi
