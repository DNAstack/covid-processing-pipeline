version 1.0

workflow assembly_to_vcf {
	input {
		String sample_id
		File assembly
		File reference_genome
		String reference_genome_id
		String date
		String snpEff_genome
		File snpEff_db
	}

	call call_variants {
		input:
			sample_id = sample_id,
			assembly = assembly,
			reference_genome = reference_genome,
			reference_genome_id = reference_genome_id
	}

	call annotate_variants {
		input:
			sample_id = sample_id,
			date = date,
			vcf = call_variants.vcf,
			reference_genome_id = reference_genome_id,
			snpEff_genome = snpEff_genome,
			snpEff_db = snpEff_db
	}

	output {
		File annotated_vcf = annotate_variants.annotated_vcf
		File annotated_vcf_index = annotate_variants.annotated_vcf_index
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
		description: "Workflow to calculate and annotate variants from SARS-CoV-2 assembly data"
	}
}

task call_variants {
	input {
		String sample_id
		File assembly
		File reference_genome
		String reference_genome_id
	}

	command {
		variants_from_assembly.sh \
			-s ~{sample_id} \
			-a ~{assembly} \
			-r ~{reference_genome} \
			-i ~{reference_genome_id}
	}

	output {
		File vcf = "~{sample_id}.vcf"
	}

	runtime {
		docker: "dnastack/variants_from_assembly:latest"
		cpu: 1
		memory: "3.75 GB"
		disks: "local-disk 20 HDD"
	}
}

task annotate_variants {
	input {
		String sample_id
		String date
		File vcf
		String reference_genome_id
		String snpEff_genome
		File snpEff_db
	}

	command {
		unzip -d /opt/snpEff/data ~{snpEff_db}

		echo -e "~{reference_genome_id}\\t1\\t29903\\t~{date}" > date.bed
		bgzip date.bed && tabix -p bed date.bed.gz

		echo -e "~{reference_genome_id}\\t~{snpEff_genome}" > to_snpeff.txt
		echo -e "~{snpEff_genome}\\t~{reference_genome_id}" > from_snpeff.txt

		vcf-annotate \
			-a date.bed.gz \
			-c "CHROM,FROM,TO,INFO/PD" \
			-d "key=INFO,ID=PD,Number=1,Type=String,Description='Processing date (YYYY-MM-DD_HH:MM:SS) (UTC)'" \
			~{vcf} | \
		bcftools annotate \
			--no-version \
			--rename-chrs to_snpeff.txt \
			- | \
		java -jar /opt/snpEff/snpEff.jar \
			ann \
			~{snpEff_genome} | \
		bcftools annotate \
			--no-version \
			--rename-chrs from_snpeff.txt \
			- \
		> ~{sample_id}.ann.vcf

		bgzip ~{sample_id}.ann.vcf
		tabix ~{sample_id}.ann.vcf.gz
	}

	output {
		File annotated_vcf = "~{sample_id}.ann.vcf.gz"
		File annotated_vcf_index = "~{sample_id}.ann.vcf.gz.tbi"
	}

	runtime {
		docker: "dnastack/variants_from_assembly:latest"
		cpu: 1
		memory: "3.75 GB"
		disks: "local-disk 20 HDD"
	}
}
