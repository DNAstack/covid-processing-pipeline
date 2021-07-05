version 1.0

workflow signal {
	input {
		String accession
		Array [File] fastq_R1s
		Array [File] fastq_R2s

		File viral_reference_genome
		File composite_reference
		String viral_reference_contig_name

		## Trimgalore, Ivar
		Int min_length
		Int min_qual

		## Ivar
		File scheme_bed
		File primer_pairs_tsv
		File amplicon_bed
		Float ivar_freq_threshold
		Int ivar_min_coverage_depth
		Float ivar_min_freq_threshold
		Int ivar_min_variant_quality

		## Samtools mpileup
		Int mpileup_depth

		## kraken
		File kraken_db

		## quast
		File viral_reference_feature_coords
	}

	call run_signal {
		input:
			accession = accession,
			fastq_R1s = fastq_R1s,
			fastq_R2s = fastq_R2s,
			min_qual = min_qual,
			min_length = min_length,
			scheme_bed = scheme_bed,
			primer_pairs_tsv = primer_pairs_tsv,
			amplicon_bed = amplicon_bed,
			composite_reference = composite_reference,
			viral_reference_genome = viral_reference_genome,
			viral_reference_feature_coords = viral_reference_feature_coords,
			viral_reference_contig_name = viral_reference_contig_name,
			kraken2_db = kraken_db,
			mpileup_depth = mpileup_depth,
			ivar_freq_threshold = ivar_freq_threshold,
			ivar_min_coverage_depth = ivar_min_coverage_depth,
			ivar_min_freq_threshold = ivar_min_freq_threshold,
			ivar_min_variant_quality = ivar_min_variant_quality
	}

	output {
		File consensus_fa = run_signal.ivar_assembly
		File vcf = run_signal.ivar_vcf
		File vcf_index = run_signal.ivar_vcf_index
		File bam = run_signal.bam
		Array [File] full_output = [run_signal.ivar_vcf, run_signal.ivar_vcf_index, run_signal.ivar_assembly, run_signal.freebayes_vcf, run_signal.freebayes_vcf_index, run_signal.freebayes_assembly, run_signal.summary, run_signal.bam]
	}
}

task run_signal {
	input {
		String accession
		Array [File] fastq_R1s
		Array [File] fastq_R2s
		Int min_qual
		Int min_length
		File scheme_bed
		File primer_pairs_tsv
		File amplicon_bed
		File composite_reference
		File viral_reference_genome
		File viral_reference_feature_coords
		String viral_reference_contig_name
		File kraken2_db
		Int mpileup_depth
		Float ivar_freq_threshold
		Int ivar_min_coverage_depth
		Float ivar_min_freq_threshold
		Int ivar_min_variant_quality
	}

	Int threads = 8
	String kraken2_db_base = basename(kraken2_db, ".tar.gz") + "/db"
	String composite_reference_base = basename(composite_reference, ".tar.gz")
	Int num_paired_fastqs = length(fastq_R1s)

	command <<<
		set -eEo pipefail

		ln -s /covid-19-signal/scripts/ "$(pwd)"

		viral_reference_contig_name='~{viral_reference_contig_name}'
		export viral_reference_contig_name

		# Generate config
		echo \
			"samples: 'sample_table.csv'
			result_dir: $(pwd)/output/
			min_qual: ~{min_qual}
			min_len: ~{min_length}
			scheme_bed: '~{scheme_bed}'
			composite_reference: '$(pwd)/~{composite_reference_base}/composite_human_viral_reference.fna'
			viral_reference_contig_name: '~{viral_reference_contig_name}'
			viral_reference_genome: '~{viral_reference_genome}'
			viral_reference_feature_coords: '~{viral_reference_feature_coords}'
			run_breseq: False
			breseq_reference: ''
			run_freebayes: False
			kraken2_db: '$(pwd)/~{kraken2_db_base}'
			primer_pairs_tsv: '-f ~{primer_pairs_tsv}'
			mpileup_depth: ~{mpileup_depth}
			var_freq_threshold: ~{ivar_freq_threshold}
			var_min_coverage_depth: ~{ivar_min_coverage_depth}
			var_min_freq_threshold: ~{ivar_min_freq_threshold}
			var_min_variant_quality: ~{ivar_min_variant_quality}
			amplicon_loc_bed: '~{amplicon_bed}'
			phylo_include_seqs: ''
			negative_control_prefix: []" \
		| tr -d '\t' > config.yaml

		tar -zxvf ~{kraken2_db}
		tar -zxvf ~{composite_reference}

		for i in {1..~{num_paired_fastqs}}; do
			echo "~{accession}" >> samplename.tmp
		done
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
			-s /covid-19-signal/Snakefile \
			all

		conda run \
			-n snakemake \
			snakemake \
			--verbose \
			--use-conda \
			--conda-prefix "$HOME/.snakemake" \
			--cores ~{threads} \
			-s /covid-19-signal/Snakefile \
			postprocess

		# iVar
		## Variants and index
		cp output/~{accession}/core/~{accession}_ivar_variants.tsv ./~{accession}.tsv
		ivar_variants_to_vcf.py ~{accession}.tsv ~{accession}.ivar.vcf
		bgzip ~{accession}.ivar.vcf && tabix ~{accession}.ivar.vcf.gz

		## Assembly
		cat output/~{accession}/core/~{accession}.consensus.fa | sed '1s/>\(.*\)/>~{accession} \1/' > ~{accession}.ivar.fa


		# Freebayes
		## Variants and index
		cp output/~{accession}/freebayes/~{accession}.variants.norm.vcf ~{accession}.freebayes.vcf
		bgzip ~{accession}.freebayes.vcf && tabix ~{accession}.freebayes.vcf.gz

		## Assembly
		cp output/~{accession}/freebayes/~{accession}.consensus.fasta ~{accession}.freebayes.fa


		# Summary info
		mkdir ~{accession}_summary
		unzip output/summary.zip -d ~{accession}_summary/
		cp ~{accession}_summary/summary.html ~{accession}_summary/~{accession}/
		cp output/*.png ~{accession}_summary/
		cp output/lineage_assignments.tsv ~{accession}_summary/
		cp output/freebayes_lineage_assignments.tsv ~{accession}_summary/
		zip -r ~{accession}_summary.zip ~{accession}_summary/
	>>>

	output {
		File ivar_vcf = "~{accession}.ivar.vcf.gz"
		File ivar_vcf_index = "~{accession}.ivar.vcf.gz.tbi"
		File ivar_assembly = "~{accession}.ivar.fa"
		File freebayes_vcf = "~{accession}.freebayes.vcf.gz"
		File freebayes_vcf_index = "~{accession}.freebayes.vcf.gz.tbi"
		File freebayes_assembly = "~{accession}.freebayes.fa"
		File summary = "~{accession}_summary.zip"
		File bam = "output/~{accession}/core/~{accession}_viral_reference.mapping.primertrimmed.sorted.bam"
	}

	runtime {
		docker: "dnastack/signal:09372f8"
		cpu: 8
		memory: "32 GB"
		disks: "local-disk 500 HDD"
		bootDiskSizeGb: 20
		preemptible: 2
	}
}
