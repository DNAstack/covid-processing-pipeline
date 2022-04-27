# COVID-19 variant calling pipeline

This repository contains workflows for processing SARS-CoV-2 data.

The [FASTQ-based workflows](#fastq-based-workflows) produce variant calls, assembled genomes, and lineage assignments from raw sequencing reads. The FASTQ-based workflows are able to process reads originating from PacBio, Oxford Nanopore (single-end), and Illumina (paired-end) sequencing data.

An [assembly-based workflow](#assembly-based-workflow) which calculates pseudo-variant sites and assigns lineage using assembled SARS-CoV-2 genomes is also included.


## Workflows

### FASTQ-based workflows

These workflows can be used to process FASTQ files into variant calls, assembled genomes, and lineage metadata.

Choose the workflow that corresponds to your sequencing data type.

The inputs and outputs for each of the FASTQ-based workflows is outlined in detail in the repository for that workflow, linked below ([PacBio](https://github.com/DNAstack/PacBio_CoSA), [Illumina](https://github.com/DNAstack/Illumina_SIGNAL), [Oxford Nanopore](https://github.com/DNAstack/Oxford_Nanopore_ARTIC)). In addition to the output files specified in those repositories, each of the FASTQ-based workflows also outputs a file containing lineage metadata calculated using [Pangolin](https://github.com/cov-lineages/pangolin) for the assembled genome that is produced during workflow execution.


#### PacBio

This workflow uses [PacBio's CoSA pipeline](https://github.com/PacificBiosciences/CoSA) to process Pacific Biosciences SARS-CoV-2 long read HiFi data.

- [More information](https://github.com/DNAstack/PacBio_CoSA)

- [Workflow inputs](https://github.com/DNAstack/PacBio_CoSA#workflow-inputs)

- [Workflow outputs](https://github.com/DNAstack/PacBio_CoSA#workflow-outputs)


#### Illumina

This workflow uses the [SIGNAL pipeline](https://github.com/jaleezyy/covid-19-signal) to process Illumina paired-end SARS-CoV-2 sequencing data.

- [More information](https://github.com/DNAstack/Illumina_SIGNAL)

- [Workflow inputs](https://github.com/DNAstack/Illumina_SIGNAL#workflow-inputs)

- [Workflow outputs](https://github.com/DNAstack/Illumina_SIGNAL#workflow-outputs)


#### Oxford Nanopore

This workflow uses the [Connor lab's implementation](https://github.com/connor-lab/ncov2019-artic-nf) of the [ARTIC pipeline](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) to process Oxford Nanopore single-ended SARS-CoV-2 sequencing data.

- [More information](https://github.com/DNAstack/Oxford_Nanopore_ARTIC)

- [Workflow inputs](https://github.com/DNAstack/Oxford_Nanopore_ARTIC#workflow-inputs)

- [Workflow outputs](https://github.com/DNAstack/Oxford_Nanopore_ARTIC#workflow-outputs)


### Assembly-based workflow

The [variants from assembly workflow](./workflows/from_assembly/assembly_to_vcf.wdl) can be used to determine 'pseudo-variant' sites when no raw sequencing data is available. This workflow aligns the provided assembled SARS-CoV-2 genome to the reference genome, then uses [`snp-sites`](https://github.com/sanger-pathogens/snp-sites) to determine sites that differ from the reference. Variant sites are output in VCF format. Viral lineage is assigned using [Pangolin](https://github.com/cov-lineages/pangolin), as in the FASTQ-based workflows.

N.B. that since base quality information is not available when using an assembly alone to call variants, these variant sites cannot be filtered based on quality and should be used for exploratory analysis only. In addition, indels cannot be called using this method. Prefer the FASTQ-based workflows when raw sequencing data is available.


#### Workflow inputs

| Input | Description |
|:-|:-|
| `accession` | Sample ID |
| `assembly` | Assembled SARS-CoV-2 genome |
| `reference_genome` | [The SARS-CoV-2 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3) |
| `reference_genome_id` | [`MN908947.3`] |
| `container_registry`  | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |


#### Workflow outputs

| Output | Description |
|:-|:-|
| `vcf`, `vcf_index` | Pseudo-variant calls and index in VCF format |
| `lineage_metadata` | Lineage assignment and associated metadata (tool versions etc.) output by `Pangolin` |


## Running workflows

### Required software

- [Docker](https://docs.docker.com/get-docker/)
- [Cromwell](https://github.com/broadinstitute/cromwell/releases) & Java (8+) OR [miniwdl](https://github.com/chanzuckerberg/miniwdl/releases) & python3

### Running using Cromwell

From the root of the repository, run:

```bash
java -jar /path/to/cromwell.jar run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in the `cromwell-executions` directory. When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.


### Running using miniwdl

This command assumes you have `miniwdl` available on your command line. If `miniwdl` is not available, try installing using `pip install miniwdl`.

```bash
miniwdl run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in a dated directory (e.g. named `20200704_073415_main`). When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.
