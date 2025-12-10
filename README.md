# Nextflow Pipeline for GATK SNP Calling

## What does this pipeline do?

This pipline works for paired end fastq.gz files that have already been downloaded. It will clean the files, filter them, map them to a reference, perform base recalibration, call germline SNPs via local re-assembly of haplotypes, create a genomic database, and perform joint genotyping. The tools used in this pipeline are gatk4 version 4.6.1.0, fastp version, 1.0.1, samtools version 1.22.1, and bwa version 0.7.19. 

## Running this pipeline using only default parameters

1) Set up your directories. In your chosen work directory (the directory you are going to run nextflow from) set up the following directories as shown by the tree command below.

``` bash
tree
.
├── BAM
├── Chr_GenomicsDB
├── data
├── Genome
|    └── Reference
├── VCF
|    └── HaplotypeCaller
├── WGS_Filtered
```

