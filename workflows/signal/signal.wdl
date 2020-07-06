version 1.0

workflow signal {
	input {
		String samplename
		Array [File] fastq_R1s
		Array [File] fastq_R2s

		File viral_reference_genome
		File human_reference

		## Trimgalore, Ivar
		Int min_length
		Int min_qual

		## Ivar
		File scheme_bed
		Float ivar_freq_threshold
		Int ivar_min_coverage_depth
		Float ivar_min_freq_threshold
		Int ivar_min_variant_quality

		Float min_freq_high_confidence_threshold

		## Samtools mpileup
		Int mpileup_depth

		## breseq
		File breseq_reference

		## kraken
		File kraken_db

		## quast
		File viral_reference_feature_coords
	}

	call run_signal {
		input:
			samplename = samplename,
			fastq_R1s = fastq_R1s,
			fastq_R2s = fastq_R2s,
			min_qual = min_qual,
			min_len = min_length,
			scheme_bed = scheme_bed,
			human_reference = human_reference,
			viral_reference_genome = viral_reference_genome,
			viral_reference_feature_coords = viral_reference_feature_coords,
			breseq_reference = breseq_reference,
			kraken2_db = kraken_db,
			mpileup_depth = mpileup_depth,
			ivar_freq_threshold = ivar_freq_threshold,
			ivar_min_coverage_depth = ivar_min_coverage_depth,
			ivar_min_freq_threshold = ivar_min_freq_threshold,
			ivar_min_variant_quality = ivar_min_variant_quality,
			min_freq_high_confidence_threshold = min_freq_high_confidence_threshold
	}

	output {
		File consensus_fa = run_signal.consensus_fa
		File vcf = run_signal.vcf
		File vcf_index = run_signal.vcf_index
		File low_freq_vcf = run_signal.low_freq_vcf
		File low_freq_vcf_index = run_signal.low_freq_vcf_index
		File bam = run_signal.bam
		File bam_index = run_signal.bam_index
		File full_output = run_signal.full_output
	}
}

task run_signal {
	input {
		String samplename
		Array [File] fastq_R1s
		Array [File] fastq_R2s
		Int min_qual
		Int min_len
		File scheme_bed
		File human_reference
		File viral_reference_genome
		File viral_reference_feature_coords
		File breseq_reference
		File kraken2_db
		Int mpileup_depth
		Float ivar_freq_threshold
		Int ivar_min_coverage_depth
		Float ivar_min_freq_threshold
		Int ivar_min_variant_quality
		Float min_freq_high_confidence_threshold
	}

	String viral_reference_base = basename(viral_reference_genome, ".fasta")
	String kraken2_db_base = basename(kraken2_db, ".tar.gz") + "/db"
	String human_reference_base = basename(human_reference, ".tar.gz")
	Int num_paired_fastqs = length(fastq_R1s)
	Int threads = 8

	command <<<
		ln -s /opt/covid-19-signal/scripts/ "$(pwd)"

		## Generate config
		echo \
		"result_dir: $(pwd)/output/
		min_qual: ~{min_qual}
		min_len: ~{min_len}
		scheme_bed: '~{scheme_bed}'
		human_reference: '$(pwd)/~{human_reference_base}'
		viral_reference_genome: '~{viral_reference_genome}'
		viral_reference_feature_coords: '~{viral_reference_feature_coords}'
		breseq_reference: '~{breseq_reference}'
		kraken2_db: '$(pwd)/~{kraken2_db_base}'
		mpileup_depth: ~{mpileup_depth}
		ivar_freq_threshold: ~{ivar_freq_threshold}
		ivar_min_coverage_depth: ~{ivar_min_coverage_depth}
		ivar_min_freq_threshold: ~{ivar_min_freq_threshold}
		ivar_min_variant_quality: ~{ivar_min_variant_quality}
		samples: 'sample_table.csv'
		phylo_include_seqs: ''" | tr -d '\t' > config.yaml

		tar -zxvf ~{kraken2_db}
		tar -zxvf ~{human_reference}

		## Generate sample_table.csv
		yes ~{samplename} | head -~{num_paired_fastqs} > samplename.tmp
		echo -e "~{sep='\n' fastq_R1s}" >> fastq_R1s.tmp
		echo -e "~{sep='\n' fastq_R2s}" >> fastq_R2s.tmp
		echo sample,r1_path,r2_path > sample_table.csv
		paste -d , samplename.tmp fastq_R1s.tmp fastq_R2s.tmp >> sample_table.csv && rm samplename.tmp fastq_R1s.tmp fastq_R2s.tmp


		conda run \
			-n snakemake \
			snakemake \
			--verbose \
			--use-conda \
			--conda-prefix "$HOME/.snakemake" \
			--cores ~{threads} \
			-s /opt/covid-19-signal/Snakefile all

		process_outputs.sh \
			-s ~{samplename} \
			-m ~{min_freq_high_confidence_threshold} \
			-v ~{viral_reference_base}

		tar -zcvf ~{samplename}.tar.gz ~{samplename}
	>>>

	output {
		File consensus_fa = "~{samplename}.consensus.fa"
		File vcf = "~{samplename}.vcf.gz"
		File vcf_index = "~{samplename}.vcf.gz.tbi"
		File low_freq_vcf = "~{samplename}.low_freq.vcf.gz"
		File low_freq_vcf_index = "~{samplename}.low_freq.vcf.gz.tbi"
		File bam = "~{samplename}.bam"
		File bam_index = "~{samplename}.bam.bai"
		File full_output = "~{samplename}.tar.gz"
	}

	runtime {
		docker: "dnastack/signal:1"
		cpu: 8
		memory: "32 GB"
		disks: "local-disk 500 HDD"
		bootDiskSizeGb: 20
		preemptible: 2
	}
}
