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
2) Download your paired end fatq.gz files into the 'data' directory. The default parameters assume that your file is named with the tags f1 and r2 for your forwrd and reverse pair reads. An example of file names would be CRR379058_f1.fq.gz and CRR379058_r2.fq.gz
3) Download main.nf and nextflow.config into your working directory. Do not change the names of these files. 
4) Make a file called chromosomes_DV10.list and put it into your work directory. This file should have all chromosome IDs in your dataset, with one ID per line. An example of what your file should look like is below.
``` bash
A01
A02
A03
A04
A05
A06
A07
A08
A09
A10
C01
C02
C03
C04
C05
C06
C07
C08
C09
```
5) In the Genome directory download your vcfs with the known sites and the corresponding index file. These files are used for base recalibration.
 ``` bash
known_sites_DarmorV10.sorted.vcf
known_sites_DarmorV10.sorted.vcf.idx
```
6) In the directory Reference download all reference files. You should have the following files:
``` bash
DarmorV10.amb
DarmorV10.ann
DarmorV10.bwt
DarmorV10_Chromosomes_Only.dict
DarmorV10_Chromosomes_Only.fa
DarmorV10_Chromosomes_Only.fa.fai
DarmorV10.dict
DarmorV10.pac
DarmorV10.sa
```
7) You can now run nextflow from your working directory.
``` bash
nextflow run main.nf
```
You can run nextflow with the additional parameters if you want better logs and traceback, which are useful if something goes wrong.
``` bash
nextflow run main.nf -with-trace -with-report -with-timeline
```
8) nextflow will copy all relevant output into the directories we set in step 1.

BAM - contains the bam files and the recalibration statistics

Chr_GenomicsDB - contains the genomic databases

VCF - contains the final VCF file

VCF/HaplotypeCaller - contains the gVCFs made by HaplotypeCaller

WGS_Filtered - contains the cleaned fastq files


nextflow will also make a directory called work, which will have all of these results as well as the logs. 

## Running this pipeline defining your own parameters
The purpose of the parameters is to make the pipline more flexible, so that the user can change the input and output file and directory names. All parameters are described below.

`reads_dir` - This describes the directory that the paired read fastq.gz files are stored in. In the above example this would be data.

`read1_tag` - This controls what nextflow expects the fastq.gz read 1 files to be called. This defualt expects files such as CRR379058_f1.fq.gz but if your read 1 files were instead named such as CRR379058_r1.fq.gz, you could set this parameter to r1 and nextflow would read your files. This pipeline will still assume that your fastq.gz file names start with the sample ID.

`read2_tag` - This controls what nextflow expects the fastq.gz read 2 files to be called. This defualt expects files such as CRR379058_r2.fq.gz but if your read 1 files were instead named such as CRR379058_read2.fq.gz, you could set this parameter to read2 and nextflow would read your files. This pipeline will still assume that your fastq.gz file names start with the sample ID.

`QC` - This parameter sets where your cleaned fatsq files are saved. In the above example this would be WGS_Filtered.

`Ref_Abbr` - This controls how multiple output files are named. The default is DV10, which represents an abbreviation of the refernce genome originally used for this pipeline. All output files would then have DV10 in their name. You can change this parameter to be an abbreviation of the reference you use or you can change it to any other identifying string you wish to use for your file names.

`BAM` - This controls where the bam files are saved. The defualt directory is BAM/

`Ref_index` - This is the path and prefix to the reference index. The defualt is /Genome/Reference/DarmorV10, these are the files described in step 6.

`reference` - This is the path to the reference fasta file and index. The defualt is /Genome/Reference/DarmorV10_Chromosomes_Only.fa, these are the other files described in step 6.

`KnownSites` - This is the path to the known sites vcf used for base recalibration. The defualt is /Genome/known_sites_DarmorV10.sorted.vcf and tehse files are described in step 5.

`chr_file` - This is the file that has a list chromosome IDs and is described in step 4. The default file name is chromosomes_DV10.list.

`VCF` - This controls where the HaplotypeCaller gVCFs are saved. The default is /VCF/HaplotypeCaller

`DB` - This controls where the genomic databases are saved. The default is /Chr_GenomicsDB

`final_VCF` - This controls where the final VCF is saved. The default is /VCF

`vcf_prefix` - This controls how the final VCF is named. The default is Samples_minicore_1-135, which would name the output VCF as Samples_minicore_1-135_whole_genome.DV10.raw.vcf.gz. This parameter changes the prefix for this file name. You can change the DV10 part of the name using the parameter Ref_Abbr.


An example of how to use the parameters when running your nextflow pipeline is below. All other parameters can be used in the same way.

``` bash
nextflow run main.nf --read1_tag r1 --read2_tag read2
```
