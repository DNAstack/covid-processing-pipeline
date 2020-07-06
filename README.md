# COVID-19 variant calling pipeline

This pipeline is used to process SARS-CoV-2 sequencing data from raw reads to variants and assembly files. This pipeline is able to process both Nanopore and Illumina paired-end sequencing data.

This workflow wraps two separate analysis pipelines: 
- The **SARS-CoV-2 Illumina GeNome Assembly Line** ([SIGNAL](https://github.com/jaleezyy/covid-19-signal)) pipeline for Illumina paired-end data; created by the [McArthur lab](http://mcarthurbioinformatics.ca/)
- An implementation of the [ARTIC protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) for use with Nanopore data by the [Connor lab](https://github.com/connor-lab/ncov2019-artic-nf)

Visit the repositories for each workflow for more information on the steps, protocols, and parameters in use.



## Required software

- [Docker](https://docs.docker.com/get-docker/)
- [Cromwell](https://github.com/broadinstitute/cromwell/releases) & Java (8+) OR [miniwdl](https://github.com/chanzuckerberg/miniwdl/releases) & python3
- cURL (for download script)


## Required data files

Some reference files are required by the SIGNAL pipeline. They can be downloaded automatically using the `download-data.sh` script in the `scripts` directory by running:

```bash
./scripts/download-data.sh
```

This will download required data files to the `data` directory at the root of the repository. These files total ~3.5 GB.



## Running the workflow

It is assumed that all commands will be run from the root of the repository directory. If run from elsewhere, the input file should be updated to use absolute file paths. Note that the `download_fastqs` task (the first task) relies on a reliable internet connection; connection to NCBI's SRA can sometimes be unstable. If the workflow fails at this step, try rerunning it at a later time when NCBI's servers may be under lighter load. 


### Selecting inputs

Inputs must be selected and filled in to a JSON-format file. A template file may be found here: [inputs/input_template.json](inputs/input_template.json).


| Input identifier                           | Description                                                          | Possible values                                                                                                                                       |
|--------------------------------------------|----------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| main.Run_ID                                | NCBI run accession ID                                        | Find SARS-CoV-2 sequencing runs [here](https://www.ncbi.nlm.nih.gov/sra/?term=txid2697049%5BOrganism:noexp%5D%20NOT%200[Mbases])                                 |
| main.type                                  | Sequencing platform used; determines which analysis pipeline is run        | Must be one of ["NANOPORE", "ILLUMINA_PE"]                                                                                                            |
| main.signal_scheme_bed                     | [SIGNAL] The BED-format primer scheme used to prepare the Illumina PE library | Only used in the ILLUMINA_PE workflow (SIGNAL). Can be the path to any properly formatted file. ARTIC primer files can be found [here](data/primer_schemes/nCoV-2019).  See section on Primer schemes below.       |
| main.artic_primer_version                  | [ARTIC] The ARTIC primer version used to prepare the Nanopore library        | Only used in the NANOPORE workflow (ARTIC). Must be one of ["V1", "V2", "V3"]. Other primer schemes for Nanopore samples are not currently supported. |
| main.human_reference                | [SIGNAL] Path to an archive file containing human reference sequences         | Downloaded using the [scripts/download-data.sh](scripts/download-data.sh) script.                                                                     |
| main.kraken_db                      | [SIGNAL] Path to an archive file containing the Kraken2 database files        | "           "            "                                                                                                                            |
| main.viral_reference_genome         | [SIGNAL] Path to the SARS-CoV-2 viral reference FASTA (MN908947.3)            | "           "            "                                                                                                                            |
| main.viral_reference_feature_coords | [SIGNAL] Path to the SARS-CoV-2 viral reference feature coordinate file       | "           "            "                                                                                                                            |
| main.breseq_reference               | [SIGNAL] Path to the SARS-CoV-2 GBK-format file                               | "           "            "                                                                                                                            |
| main.min_freq_high_confidence_threshold | [SIGNAL] Frequency threshold to call high-confidence variants             | Default: 0.25 |
| main.ivar_min_freq_threshold        | [SIGNAL] Minimum frequency threshold to call variants                         | Default: 0.03 |
| main.ivar_min_coverage_depth        | [SIGNAL] Minimum coverage depth to call variants                              | Default: 10 |
| main.ivar_min_variant_quality       | [SIGNAL] Minimum variant quality                                              | Default: 20 |
| main.ivar_freq_threshold            | [SIGNAL] Frequency threshold to build consensus sequence                      | Default: 0.75 |
| main.mpileup_depth                  | [SIGNAL] Mpileup depth for iVar                                               | Default: 100000 |
| main.signal_min_qual                       | [SIGNAL] Minimum quality used in trimgalore and iVar trim                     | Default: 20 |
| main.signal_min_length                     | [SIGNAL] Minimum read length to retain after trimming                         | Default: 20 |
| main.artic_min_length                      | [ARTIC] Minimum read length for artic guppyplex                              | Default: 400 |
| main.artic_max_length                      | [ARTIC] Maximum read length for artic guppyplex                              | Default: 700 |


All [SIGNAL] file locations have been prefilled in the [template input file](inputs/input_template.json). These files are downloaded and moved to the correct locations by the data download script.

`Run_ID`, `type`, `signal_scheme_bed` file location, and `artic_primer_version` must be filled in on the template input file for the workflow to be run. For formatting help, check out the test input files in the same directory. Any inputs with default values defined to not need to be explicitly defined in the input file if no change to the default value is required.

For more information about SIGNAL parameters, see [here](https://github.com/jaleezyy/covid-19-signal/blob/master/example_config.yaml). For more information about ARTIC parameters, see [here](https://github.com/connor-lab/ncov2019-artic-nf/blob/master/modules/help.nf).


#### Test input files

- Example Illumina PE input file: [inputs/illumina_pe.test_input.json](inputs/illumina_pe.test_input.json)
- Example Nanopore input file: [inputs/nanopore.test_input.json](inputs/nanopore.test_input.json)


#### Primer schemes

Primer schemes (the input to `main.signal_scheme_bed`) available in the [data/primer_schemes](data/primer_schemes) directory are from the official [artic-network github](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019). The v3 scheme was altered such that it's 5th column was set to 60 for every row to conform with the BED specification (required by ivar).

Different schemes may be used for the Illumina PE workflow, e.g. [schemes available in the SIGNAL github](https://github.com/jaleezyy/covid-19-signal/tree/master/resources/primer_schemes). The Nanopore workflow is currently only able to operate on ARTIC primer schemes (v1, v2, and v3).




### Running using Cromwell

*Ensure you have first filled out all input files appropriately.*

From the root of the repository, run:

```bash
java -jar /path/to/cromwell.jar run workflows/main.wdl -i inputs/input_template.json
```

Output and execution files will be located in the `cromwell-executions` directory. When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.



### Running using miniwdl

*Ensure you have first filled out all input files appropriately.*

This command assumes you have `miniwdl` available on your command line. If `miniwdl` is not available, try installing using `pip install miniwdl`.

```bash
miniwdl run workflows/main.wdl -i inputs/input_template.json
```

Output and execution files will be located in a dated directory (e.g. named `20200704_073415_main`). When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.



## Outputs

| Output identifier                            | Description                                                                        | Notes                                                         |
|----------------------------------------------|------------------------------------------------------------------------------------|---------------------------------------------------------------|
| raw_reads                                    | Fastq file(s)                                                                      | 1 file for Nanopore runs; 2 for Illumina paired-end runs.     |
| consensus_fa                                 | Assembled genome consensus sequence                                                |                                                               |
| vcf, vcf_index                               | High confidence variant calls and index                                            |                                                               |
| low_confidence_vcf, low_confidence_vcf_index | Low confidence (but more numerous) variant calls and index                         |                                                               |
| bam, bam_index                               | Read alignment to the SARS-CoV-2 reference and index                               |                                                               |
| full_output                                  | Full output from each pipeline; differs between Nanopore and Illumina_PE pipelines | Contains additional quality metrics, trimmed alignments, etc. |
