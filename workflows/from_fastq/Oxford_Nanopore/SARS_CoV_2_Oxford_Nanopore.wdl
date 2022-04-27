version 1.0

import "https://raw.githubusercontent.com/DNAstack/Oxford_Nanopore_ARTIC/main/workflows/artic.wdl" as artic
import "https://raw.githubusercontent.com/DNAstack/covid-processing-pipeline/master/workflows/common/common.wdl" as common

workflow SARS_CoV_2_Oxford_Nanopore {
	input {
		String accession
		Array[File] fastqs
		String primer_version
		Int min_length
		Int max_length
		Int min_reads

		String container_registry
	}

	call artic.run_artic {
		input:
			accession = accession,
			fastqs = fastqs,
			primer_version = primer_version,
			min_length = min_length,
			max_length = max_length,
			min_reads = min_reads,
			container_registry = container_registry
	}

	call common.assign_lineage {
		input:
			accession = accession,
			assembly = run_artic.consensus_fa,
			container_registry = container_registry
	}

	output {
		File consensus_fa = run_artic.consensus_fa
		File vcf = run_artic.vcf
		File vcf_index = run_artic.vcf_index
		File bam = run_artic.bam
		File summary = run_artic.summary
		File lineage_metadata = assign_lineage.lineage_metadata
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}
