version 1.0

import "signal/signal.wdl" as illumina_paired
import "artic/artic.wdl" as nanopore

workflow main {
	input {
		String Run_ID
		String type

		# ARTIC
		String artic_primer_version
		Int artic_min_length = 400
		Int artic_max_length = 700

		# SIGNAL
		File signal_scheme_bed

		File viral_reference_genome
		File human_reference

		## Trimgalore, Ivar
		Int signal_min_length = 20
		Int signal_min_qual = 20

		## Ivar
		Float ivar_freq_threshold = 0.75
		Int ivar_min_coverage_depth = 10
		Float ivar_min_freq_threshold = 0.03
		Int ivar_min_variant_quality = 20

		Float min_freq_high_confidence_threshold = 0.25

		## Samtools mpileup
		Int mpileup_depth = 100000

		## breseq
		File breseq_reference

		## kraken
		File kraken_db

		## quast
		File viral_reference_feature_coords
	}

	call download_fastqs {
		input:
			Run_ID = Run_ID,
			type = type
	}

	if (type == "ILLUMINA_PE") {
		call illumina_paired.signal {
			input:
				samplename = Run_ID,
				fastq_R1s = download_fastqs.R1,
				fastq_R2s = download_fastqs.R2,
				scheme_bed = signal_scheme_bed,
				min_qual = signal_min_qual,
				min_length = signal_min_length,
				human_reference = human_reference,
				viral_reference_genome = viral_reference_genome,
				viral_reference_feature_coords = viral_reference_feature_coords,
				breseq_reference = breseq_reference,
				kraken_db = kraken_db,
				mpileup_depth = mpileup_depth,
				ivar_freq_threshold = ivar_freq_threshold,
				ivar_min_coverage_depth = ivar_min_coverage_depth,
				ivar_min_freq_threshold = ivar_min_freq_threshold,
				ivar_min_variant_quality = ivar_min_variant_quality,
				min_freq_high_confidence_threshold = min_freq_high_confidence_threshold	
		}
	}

	if (type == "NANOPORE") {
		call nanopore.artic {
			input:
				samplename = Run_ID,
				fastqs = download_fastqs.SE_read,
				artic_primer_version = artic_primer_version,
				min_length = artic_min_length,
				max_length = artic_max_length
		}
	}

	output {
		Array [File] raw_reads = flatten([download_fastqs.R1, download_fastqs.R2, download_fastqs.SE_read_zipped])
		File consensus_fa = select_first([signal.consensus_fa, artic.consensus_fa])
		File vcf = select_first([signal.vcf, artic.vcf])
		File vcf_index = select_first([signal.vcf_index, artic.vcf_index])
		File low_confidence_vcf = select_first([signal.low_freq_vcf, artic.unfiltered_vcf])
		File low_confidence_vcf_index = select_first([signal.low_freq_vcf_index, artic.unfiltered_vcf_index])
		File bam = select_first([signal.bam, artic.bam])
		File bam_index = select_first([signal.bam_index, artic.bam_index])
		File full_output = select_first([signal.full_output, artic.full_output])
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
		description: "Workflows for variant calling in Nanopore and Illumina paired-end SARS-CoV-2 sequencing data"
	}
}

task download_fastqs {
	input {
		String Run_ID
		String type
	}

	Int threads = 8

	command {
		retries=5
		status=1
		reads_invalid=""
		while [ $retries -gt 0 ] && [ $status -ne 0 ] || [ -n "$reads_invalid" ]; do
			fasterq-dump \
				--threads ~{threads} \
				--mem 6GB \
				--split-3 \
				~{Run_ID} \
				2> dump.stderr.txt
			status=$?
			reads_invalid=$(grep "reads invalid" dump.stderr.txt)
			((retries--))
		done

		if [ $status -ne 0 ] || [ -n "$reads_invalid" ]; then
			>&2 echo -e "[ERROR]: Error downloading fastqs.\nrc: [$status]\ninvalid reads: [$reads_invalid]"
			exit 1
		fi

		if [ ~{type} = "ILLUMINA_PE" ]; then
			if [ -e ~{Run_ID}_1.fastq ]; then
				gzip ~{Run_ID}_1.fastq
			else
				echo "ERROR: Could not find read 1 for ILLUMINA_PE sample"
				exit 1
			fi
	
			if [ -e ~{Run_ID}_2.fastq ]; then
				gzip ~{Run_ID}_2.fastq
			else
				echo "ERROR: Could not find read 2 for ILLUMINA_PE sample"
				exit 1
			fi
		elif [ ~{type} = "NANOPORE" ]; then
			if [ -e ~{Run_ID}.fastq ]; then
				gzip --keep ~{Run_ID}.fastq
			else
				echo "ERROR: Could not find SE reads for NANOPORE sample"
				exit 1
			fi
		fi
	}

	output {
		Array [File] R1 = glob("*_1.fastq.gz")
		Array [File] R2 = glob("*_2.fastq.gz")
		Array [File] SE_read = glob("~{Run_ID}.fastq")
		Array [File] SE_read_zipped = glob("~{Run_ID}.fastq.gz")
	}

	runtime {
		docker: "dnastack/sra-toolkit:2.10.7"
		cpu: threads
		memory: "7.20 GB"
		disks: "local-disk 50 HDD"
	}
}
