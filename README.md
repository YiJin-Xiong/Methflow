# Methflow

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**methflow** is a bioinformatics analysis pipeline for Methylation (Nanopore DNA/RNA) sequencing data that can be used to perform basecalling, alignment, modcalling, SNV calling, phasing and DMR analysis. 



![methflow(v1.0.0)](docs/images/methflow(v1.0.0).jpg)



The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible.  The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a [full-sized dataset](https://github.com/nf-core/test-datasets/tree/nanoseq#full-sized-test-data) obtained from the [Singapore Nanopore Expression Consortium](https://github.com/GoekeLab/sg-nex-data) on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/nanoseq/results).

## Pipeline Summary

1. Basecall ([`dorado`](https://github.com/nanoporetech/dorado))
   - Call nucleotide bases (ACGT) from POD5 files
2. Alignment ([`minimap2`](https://github.com/lh3/minimap2) (default) or [`Winnowmap`](https://github.com/marbl/Winnowmap))
   - Align genomic reads to reference genome
   - minimap2/Winnowmap workflow:
     1. Convert BAM to FASTQ
     2. Alignment
     3. Convert SAM to BAM
     4. Sort
3. Modcall ([`DeepSignal3`](https://github.com/PengNi/deepsignal3) (default) and/or [`rockfish`](https://github.com/lbcb-sci/rockfish))
   - Detect DNA methylation
4. SNV calling ([`Clair3`](https://github.com/HKU-BAL/Clair3))
   - Call germline small variant, such as SNP, insertions or deletions
5. Phase ([`WhatsHap`](https://github.com/whatshap/whatshap))
   - Phase genomic variant
6. DMR ([`Modkit`](https://github.com/nanoporetech/modkit) (default) or [`DSS`](https://bioconductor.org/packages/release/bioc/html/DSS.html) or [`PoreMeth2`](https://github.com/Lab-CoMBINE/PoreMeth2))
   - Differential methylation analysis


## Installation

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`) .
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)) for full pipeline reproducibility. 
3. Install [`Conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) .

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.tsv`:

```tsv
Group_ID	Sample_ID  	Analyte_Type	Pore_Type	Path	Genome	Species	Sample_Rate
hg002	hg002_subset1	DNA	r10	./demo/hg002/subset_1.pod5	./demo/chm13v2.0.fa	Human	5k
hg002	hg002_subset2	DNA	r10	./demo/hg002/subset_2.pod5	./demo/chm13v2.0.fa	Human	5k
```

> Each row represents a pod5 file.
>
> - **Group_ID**: Group identifier of the sample. This will be identical for replicate samples from the same experimental group.
> - **Sample_ID**: Name of the sample.
> - **Analyte_Type**: Type of the sample, should be DNA or RNA.
> - **Pore_Type**: Nanopore type of the sample, should be r9 or r10.
> - **Path**: Path of POD5 file. 
> - **Genome**: Path of FASTA genome file.
> - **Species**: Species of the sample, should be human or plant.
> - **Sample_Rate**: Sample rate of the sample during sequencing, should be 4k or 5k.

An example input samplesheet can be found [here](samplesheet.tsv).

Now, you can run the pipeline using default parameters as:

```bash
nextflow run main.nf \
    --input samplesheet.tsv \
    --protocol DNA \
    --profile <docker/singularity>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation]( ) and the [parameter documentation]( ).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results]( ) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the [output documentation]( ).

Folder structures of the outputs:

```text

├── results/sampleName
    │
    ├── basecall
    │   └── basecall.bam
    │
    ├── alignment
    │   ├── basecall.fastq
    │   ├── aligned.bam
    │   ├── aligned.sam
    │   └── aligned_sorted.bam
    │
    ├── modcall
    │   ├── pod5.CG.call_mods.tsv (deepsignal3)
    │   ├── pod5.CG.call_mods.frequency.tsv (deepsignal3)
    │   ├── predictions.tsv (rockfish)
    │   └── 
    │
    ├── SNVcall
    │   ├── clair3_merge.vcf.gz
    │   └── clair3_pileup.vcf.gz.tbi
    │
    ├── phase
    │   ├── SNV_PASS.vcf
    │   ├── phased.vcf.gz
    │   ├── phased.vcf.gz.tbi
    │   ├── whatshap_haplotagged.bam
    │   ├── whatshap_haplotagged.bam.bai
    │   ├── whatshap_haplotagged.readlist
    │   └── 
    │
    ├── DMR
    │   ├── haplotype
    │   │   ├── *.bed
    │   │   │
    │   ├── PoreMeth2/ 
    │   │   ├── extracted_sorted.entropy.file.tsv
    │   │   ├── ExpQualityPlot.png
    │   │   ├── DMRTable.txt
    │   │   ├── DMRAnnotatedTable.txt
    │   │   └──  poremeth2_log.log
    │   │   │
    │   │   ├── DSS/ 
    │   │   ├── haplotype_preprocessed_1.bed
    │   │   ├── haplotype_preprocessed_2.bed
    │   │   ├── DMLtest.txt
    │   │   ├── callDML.txt
    │   │   ├── callDMR.txt
    │   │   └── logdmr.log
    │   │
    │   ├── modkit_dmr.bed
    │   ├── modkit_dmr.log
        └── 

```

- pipeline_info: Information of the workflow execution
- sampleName: directory to save the results of Methflow, name set by `Sample_ID` in `samplesheet.tsv`
  - basecall: basecall files
  - alignment: alignment files
  - modcall: modcall results
  - SNVcall: SNVs generated by Clair3
  - phase: Phased SNVs by WhatsHap
  - DMR: DMLs and DMRs generated by modkit or PoreMeth2



## Acknolwedgements

- Some code were referenced from [nanome](https://github.com/TheJacksonLaboratory/nanome) and [nf-core](https://github.com/nf-core).
- Code for ASM detection using DSS were referenced from [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase) of Akbari *et al.*
