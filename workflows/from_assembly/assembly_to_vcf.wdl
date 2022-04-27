version 1.0

import "https://raw.githubusercontent.com/DNAstack/covid-processing-pipeline/master/workflows/common/common.wdl" as common

workflow assembly_to_vcf {
	input {
		String accession
		File assembly
		File reference_genome
		String reference_genome_id

		String container_registry
	}

	call call_variants {
		input:
			accession = accession,
			assembly = assembly,
			reference_genome = reference_genome,
			reference_genome_id = reference_genome_id,
			container_registry = container_registry
	}

	call common.assign_lineage {
		input:
			accession = accession,
			assembly = assembly,
			container_registry = container_registry
	}

	output {
		File vcf = call_variants.vcf
		File vcf_index = call_variants.vcf_index
		File lineage_metadata = assign_lineage.lineage_metadata
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}

task call_variants {
	input {
		String accession
		File assembly
		File reference_genome
		String reference_genome_id

		String container_registry
	}

	command {
		variants_from_assembly.sh \
			-s ~{accession} \
			-a ~{assembly} \
			-r ~{reference_genome} \
			-i ~{reference_genome_id}

		bgzip "~{accession}.vcf"
		tabix "~{accession}.vcf.gz"
	}

	output {
		File vcf = "~{accession}.vcf.gz"
		File vcf_index = "~{accession}.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/variants_from_assembly:latest"
		cpu: 1
		memory: "3.75 GB"
		disks: "local-disk 20 HDD"
	}
}
