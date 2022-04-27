version 1.0

import "https://raw.githubusercontent.com/DNAstack/Illumina_SIGNAL/main/workflows/signal.wdl" as signal
import "https://raw.githubusercontent.com/DNAstack/covid-processing-pipeline/master/workflows/common/common.wdl" as common

workflow SARS_CoV_2_Illumina {
	input {
		String accession
		Array[File] fastq_R1s
		Array[File] fastq_R2s
		File scheme_bed
		File viral_reference_genome
		File viral_reference_feature_coords
		String viral_reference_contig_name
		File primer_pairs_tsv
		File amplicon_bed

		String container_registry
	}

	call signal.run_signal {
		input:
			accession = accession,
			fastq_R1s = fastq_R1s,
			fastq_R2s = fastq_R2s,
			scheme_bed = scheme_bed,
			viral_reference_genome = viral_reference_genome,
			viral_reference_feature_coords = viral_reference_feature_coords,
			viral_reference_contig_name = viral_reference_contig_name,
			primer_pairs_tsv = primer_pairs_tsv,
			amplicon_bed = amplicon_bed,
			container_registry = container_registry
	}

	call common.assign_lineage {
		input:
			accession = accession,
			assembly = run_signal.ivar_assembly,
			container_registry = container_registry
	}

	output {
		File ivar_vcf = run_signal.ivar_vcf
		File ivar_vcf_index = run_signal.ivar_vcf_index
		File ivar_assembly = run_signal.ivar_assembly
		File freebayes_vcf = run_signal.freebayes_vcf
		File freebayes_vcf_index = run_signal.freebayes_vcf_index
		File freebayes_assembly = run_signal.freebayes_assembly
		File summary = run_signal.summary
		File signal_lineage_metadata = run_signal.lineage_metadata
		File bam = run_signal.bam
		File lineage_metadata = assign_lineage.lineage_metadata
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}
