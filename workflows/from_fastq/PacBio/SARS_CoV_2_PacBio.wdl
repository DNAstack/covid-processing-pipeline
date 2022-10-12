version 1.0

import "https://raw.githubusercontent.com/DNAstack/PacBio/main/CoSA-CoronavirusSequencingAnalysis/workflows/cosa.wdl" as cosa
import "https://raw.githubusercontent.com/DNAstack/covid-processing-pipeline/master/workflows/common/common.wdl" as common

workflow SARS_CoV_2_PacBio {
	input {
		String accession
		File primer_trimmed_fastq

		File reference
		File reference_index

		String container_registry
	}

	call cosa.run_cosa {
		input:
			accession = accession,
			primer_trimmed_fastq = primer_trimmed_fastq,
			reference = reference,
			reference_index = reference_index,
			container_registry = container_registry
	}

	call common.assign_lineage {
		input:
			accession = accession,
			assembly = run_cosa.consensus_fa,
			container_registry = container_registry
	}

	output {
		File aligned_bam = run_cosa.aligned_bam
		File aligned_bam_index = run_cosa.aligned_bam_index
		File high_quality_variants = run_cosa.high_quality_variants
		File high_quality_variants_index = run_cosa.high_quality_variants_index
		File consensus_fa = run_cosa.consensus_fa
		File lineage_metadata = assign_lineage.lineage_metadata
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}
