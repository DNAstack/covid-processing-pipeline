version 1.0

workflow artic {
	input {
		String samplename
		Array [File] fastqs
		String artic_primer_version
		Int min_length
		Int max_length
	}

	call run_artic {
		input:
			samplename = samplename,
			fastqs = fastqs,
			artic_primer_version = artic_primer_version,
			min_length = min_length,
			max_length = max_length
	}

	output {
		File vcf = run_artic.vcf
		File vcf_index = run_artic.vcf_index
		File consensus_fa = run_artic.consensus_fa
		File unfiltered_vcf = run_artic.unfiltered_vcf
		File unfiltered_vcf_index = run_artic.unfiltered_vcf_index
		File bam = run_artic.bam
		File bam_index = run_artic.bam_index
		File full_output = run_artic.full_output
	}
}

task run_artic {
	input {
		String samplename
		Array [File] fastqs
		String artic_primer_version
		Int min_length
		Int max_length
	}

	command <<<
		cp ~{sep=' ' fastqs} .

		nextflow run /opt/ncov2019-artic-nf/main.nf \
			-profile conda \
			--cache /opt/conda/envs/ \
			--medaka \
			--prefix ~{samplename} \
			--basecalled_fastq . \
			--schemeRepoURL /opt/artic-ncov2019 \
			--schemeVersion ~{artic_primer_version} \
			--min_length ~{min_length} \
			--max_length ~{max_length} \
			--outdir ~{samplename}

		qc_pass=$(tail -1 ~{samplename}/~{samplename}.qc.csv | awk -F ',' '{print $NF}')
		if [ ! "$qc_pass" = "TRUE" ]; then
			exit 1
		fi

		tar -zcvf ~{samplename}.tar.gz ~{samplename}

		process_outputs.sh \
			-s ~{samplename}
	>>>

	output {
		File consensus_fa = "~{samplename}.consensus.fa"
		File vcf = "~{samplename}.vcf.gz"
		File vcf_index = "~{samplename}.vcf.gz.tbi"
		File unfiltered_vcf = "~{samplename}.medaka.vcf.gz"
		File unfiltered_vcf_index = "~{samplename}.medaka.vcf.gz.tbi"
		File bam = "~{samplename}.bam"
		File bam_index = "~{samplename}.bam.bai"
		File full_output = "~{samplename}.tar.gz"
	}

	runtime {
		docker: "dnastack/artic:1"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 200 HDD"
		preemptible: 2
	}
}