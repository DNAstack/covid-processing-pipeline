version 1.0

import "signal/signal.wdl" as illumina_paired
import "artic/artic.wdl" as nanopore

workflow main {
	input {
		String accession
		String type

		# ARTIC
		String artic_primer_version
		Int artic_min_length = 400
		Int artic_max_length = 700

		# SIGNAL
		File signal_scheme_bed
		File signal_primer_pairs_tsv
		File signal_amplicon_bed

		File viral_reference_genome
		File composite_reference
		String viral_reference_contig_name

		## Trimgalore, Ivar
		Int signal_min_length = 20
		Int signal_min_qual = 20

		## Ivar
		Float ivar_freq_threshold = 0.75
		Int ivar_min_coverage_depth = 10
		Float ivar_min_freq_threshold = 0.03
		Int ivar_min_variant_quality = 20

		## Samtools mpileup
		Int mpileup_depth = 100000

		## kraken
		File kraken_db

		## quast
		File viral_reference_feature_coords
	}

	call download_fastqs {
		input:
			accession = accession,
			type = type
	}

	if (type == "ILLUMINA_PE") {
		call illumina_paired.signal {
			input:
				accession = accession,
				fastq_R1s = download_fastqs.R1,
				fastq_R2s = download_fastqs.R2,
				scheme_bed = signal_scheme_bed,
				primer_pairs_tsv = signal_primer_pairs_tsv,
				amplicon_bed = signal_amplicon_bed,
				min_qual = signal_min_qual,
				min_length = signal_min_length,
				composite_reference = composite_reference,
				viral_reference_genome = viral_reference_genome,
				viral_reference_feature_coords = viral_reference_feature_coords,
				viral_reference_contig_name = viral_reference_contig_name,
				kraken_db = kraken_db,
				mpileup_depth = mpileup_depth,
				ivar_freq_threshold = ivar_freq_threshold,
				ivar_min_coverage_depth = ivar_min_coverage_depth,
				ivar_min_freq_threshold = ivar_min_freq_threshold,
				ivar_min_variant_quality = ivar_min_variant_quality
		}
	}

	if (type == "NANOPORE") {
		call nanopore.artic {
			input:
				samplename = accession,
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
		File bam = select_first([signal.bam, artic.bam])
		Array [File] full_output = select_first([signal.full_output, artic.full_output])
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
		description: "Workflows for variant calling in Nanopore and Illumina paired-end SARS-CoV-2 sequencing data"
	}
}

task download_fastqs {
	input {
		String accession
		String type
	}

	command {
		set -Eo pipefail

		use_fastq_dump=false

		retries=5
		status=1

		fastq_dump=""
		reads_invalid=""

		while [ $retries -gt 0 ] && [ $status -ne 0 ] || [ -n "$reads_invalid" ]; do
			echo "Trying download with [$retries] retries (fasterq-dump)"

			fasterq-dump \
				--mem 3GB \
				--split-3 \
				~{accession} \
				2> dump.stderr.txt
			status=$?
			reads_invalid=$(grep "reads invalid" dump.stderr.txt)
			fastq_dump=$(grep "fastq-dump" dump.stderr.txt)
			if [ -n "$fastq_dump" ]; then
				cat dump.stderr.txt
				use_fastq_dump=true
				break
			fi
			((retries--))
		done

		# Try using fastq-dump rather than fasterq-dump
		if [ "$use_fastq_dump" = "true" ] || [ $status -ne 0 ] || [ -n "$reads_invalid" ]; then
			retries=5
			status=1
			reads_invalid=""
			while [ $retries -gt 0 ] && [ $status -ne 0 ] || [ -n "$reads_invalid" ]; do
				echo "Trying download with [$retries] retries (fastq-dump)"

				fastq-dump \
					--split-e \
					~{accession} \
					2> dump.stderr.txt \
				status=$?
				reads_invalid=$(grep "reads invalid" dump.stderr.txt)
				((retries--))
			done
		fi

		if [ $status -ne 0 ] || [ -n "$reads_invalid" ]; then
			>&2 echo -e "[ERROR]: Error downloading fastqs.\\nrc: [$status]\\ninvalid reads: [$reads_invalid]\\nStderr:\\n$(cat dump.stderr.txt)"
			exit 1
		fi

		if [ ~{type} = "ILLUMINA_PE" ]; then
			if [ -e ~{accession}_1.fastq ]; then
				bgzip ~{accession}_1.fastq
			else
				echo "ERROR: Could not find read 1 for ILLUMINA_PE sample"
				exit 1
			fi

			if [ -e ~{accession}_2.fastq ]; then
				bgzip ~{accession}_2.fastq
			else
				echo "ERROR: Could not find read 2 for ILLUMINA_PE sample"
				exit 1
			fi
		elif [ ~{type} = "NANOPORE" ]; then
			if [ -e ~{accession}.fastq ]; then
				bgzip -c ~{accession}.fastq > ~{accession}.fastq.gz
			else
				echo "ERROR: Could not find SE reads for NANOPORE sample"
				exit 1
			fi
		fi
	}

	output {
		Array [File] R1 = glob("*_1.fastq.gz")
		Array [File] R2 = glob("*_2.fastq.gz")
		Array [File] SE_read = glob("~{accession}.fastq")
		Array [File] SE_read_zipped = glob("~{accession}.fastq.gz")
	}

	runtime {
		docker: "dnastack/sra-toolkit-tabix:2.10.7"
		cpu: 1
		memory: "3.5 GB"
		disks: "local-disk 50 HDD"
	}
}
