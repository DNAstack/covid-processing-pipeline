version 1.0

task assign_lineage {
	input {
		String accession
		File assembly

		String container_registry
	}

	Int threads = 1

	command {
		pangolin \
			--threads ~{threads} \
			--outfile ~{accession}.lineage_info.csv \
			~{assembly}

		assembly_length=$(sed '1d' ~{assembly} | tr -d '\n' | sed -e 's/^N*//' -e 's/N*$//' | wc -c)
		paste -d ',' ~{accession}.lineage_info.csv <(echo -e "assembly_length\n$assembly_length") > ~{accession}.lineage.csv
	}

	output {
		File lineage_metadata = "~{accession}.lineage.csv"
	}

	runtime {
		docker: "~{container_registry}/pangolin:20eb73e"
		cpu: threads
		memory: "7.5 GB"
		disks: "local-disk 50 HDD"
	}
}
