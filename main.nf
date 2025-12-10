params.reads_dir = "data"
params.read1_tag = "f1"
params.read2_tag = "r2"
params.QC = "${baseDir}/WGS_Filtered"
params.Ref_Abbr = "DV10"
params.BAM = "${baseDir}/BAM"
params.Ref_index = "${baseDir}/Genome/Reference/DarmorV10"
params.reference   = "${baseDir}/Genome/Reference/DarmorV10_Chromosomes_Only.fa"
params.KnownSites = "${baseDir}/Genome/known_sites_DarmorV10.sorted.vcf"
params.chr_file = "${baseDir}/chromosomes_DV10.list"
params.VCF = "${baseDir}/VCF/HaplotypeCaller"
params.DB = "${baseDir}/Chr_GenomicsDB"
params.final_VCF = "${baseDir}/VCF"
params.vcf_prefix = "Samples_minicore_1-135"

process fastp {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.QC}", mode: 'copy'

    input:
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id),
          path("${read1.baseName}_clean.fastq.gz"),
          path("${read2.baseName}_clean.fastq.gz"),
          path("fastp_${id}_*.{html,json}")

    script:
    """
    OUT1=${read1.baseName}_clean.fastq.gz
    OUT2=${read2.baseName}_clean.fastq.gz

    ACCESSION=$id

    fastp -i "$read1" -I "$read2" -o "\$OUT1" -O "\$OUT2" -h "fastp_\${ACCESSION}_${task.index}.html" -j "fastp_\${ACCESSION}_${task.index}.json" --correction --dedup --dup_calc_accuracy 5

    echo ""
    echo "Processed \$ACCESSION"
    echo ""
    echo "In: $read1"
    echo "In: $read2"
    echo ""
    echo "Out1: \$OUT1"
    echo "Out2: \$OUT2"
    echo ""
    """
}

process mapping {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.BAM}", mode: 'copy'

    input:
    tuple val(id), path(clean_read1), path(clean_read2)
    path ref_index_files

    output: tuple val(id), 
          path("${id}.${params.Ref_Abbr}.sorted.bam")

    script:
    """
    ACCESSION=$id
    SAM=\${ACCESSION}.sam

    echo ""
    echo "Processing sample --> \$ACCESSION"
    echo "Read1: $clean_read1"
    echo "Read2: $clean_read2"
    echo ""
    echo "Starting mapping"
    echo ""

    bwa mem -M -t 64 -R "@RG\\tID:\${ACCESSION}\\tLB:LIB-1\\tPL:ILLUMINA\\tSM:\${ACCESSION}\\tPU:\${ACCESSION}" "${params.Ref_index}" "$clean_read1" "$clean_read2" > "\$SAM"

    echo ""
    echo "Sorting and Indexing ..."
    echo ""

    samtools sort -@ 64 "\$SAM" -o "\${ACCESSION}.${params.Ref_Abbr}.sorted.bam"

    echo ""
    echo "Done!"
    echo ""
    """
}

process DupValIndx {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.BAM}", mode: 'copy'

    input:
    tuple val(id), path(sorted_bam)

    output: 
    tuple val(id),
          path("${id}.${params.Ref_Abbr}.marked.bam"),
          path("${id}.${params.Ref_Abbr}.marked*.bai"),
          path("${id}.${params.Ref_Abbr}.metrics.txt"),
          path("${id}.${params.Ref_Abbr}.validation.txt")

    script:
    """
    DEDUP_BAM="${id}.${params.Ref_Abbr}.marked.bam"
    METRICS="${id}.${params.Ref_Abbr}.metrics.txt"
    VALIDATION="${id}.${params.Ref_Abbr}.validation.txt"
    METRICS_PREFIX="${id}_${params.Ref_Abbr}_mapping_metrics"

    echo ""
    echo "Processing sample --> ${id}"
    echo "Input BAMs (coordinate sorted): ${sorted_bam}"
    echo ""
    echo "Marking Duplicates ..."
    echo ""

    gatk --java-options "-Xmx100g" MarkDuplicates -I "${sorted_bam}" -O "\$DEDUP_BAM" -M "\$METRICS"

    echo ""
    echo "Done!"
    echo ""
    echo "Validating BAM file"
    echo ""

    gatk ValidateSamFile -I "\$DEDUP_BAM" -O "\$VALIDATION" --MODE SUMMARY

    echo ""
    echo "Done!"
    echo ""
    echo "Build BAM index"
    echo ""

    gatk BuildBamIndex -I "\$DEDUP_BAM" 

    echo ""
    echo "Done!"
    echo ""
    """
    }

process recalibrate {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.BAM}", mode: 'copy'

    input:
    tuple val(id), path(marked_bam), path(reference_dir), path(known_sites)

    output: 
    tuple val(id),
          path("${id}.${params.Ref_Abbr}.recal.marked.bam"),
          path("${id}.${params.Ref_Abbr}.recal.marked*bai"),
          path("${id}_${params.Ref_Abbr}_Recal_Data.table"),
          path("${id}.${params.Ref_Abbr}.recal.validation.txt")

    script:
    """
    RECAL_TABLE="${id}_${params.Ref_Abbr}_Recal_Data.table"
    RECAL_BAM="${id}.${params.Ref_Abbr}.recal.marked.bam"
    VALIDATION="${id}.${params.Ref_Abbr}.recal.validation.txt"


    echo "Started at: \$(date)"
    echo ""
    echo "Processing sample --> ${id}"
    echo "Input BAM: ${marked_bam}"
    echo ""
    echo "Recalibrating Base Score Quality ..."
    echo ""

    gatk BaseRecalibrator -R "${params.reference}" -I "${marked_bam}" --known-sites "${params.KnownSites}" -O "\$RECAL_TABLE"

    echo ""
    echo "Done!"
    echo ""
    echo "Applying BQSR ..."
    echo ""

    gatk ApplyBQSR -R "${params.reference}" -I "${marked_bam}" --bqsr-recal-file "\$RECAL_TABLE" -O "\$RECAL_BAM"

    echo ""
    echo "Done!"
    echo ""
    echo "Validating BAM file"
    echo ""

    gatk ValidateSamFile -I "\$RECAL_BAM" -O "\$VALIDATION" --MODE SUMMARY

    echo ""
    echo "Done!"
    echo ""
    echo "Build BAM index"
    echo ""

    gatk BuildBamIndex -I "\$RECAL_BAM"

    echo ""
    echo "Done!"
    echo ""
    echo "Finished at: \$(date)"
    echo ""
    """
    }

process BQSR {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.BAM}", mode: 'copy'

    input:
    tuple val(id), path(recal_bam), path(data_table), path(reference_dir), path(known_sites)

    output: 
    tuple val(id),
          path("${id}_${params.Ref_Abbr}_Recal_Data.after.table"),
          path("${id}_Recalibration.csv")


    script:
    """
    echo ""
    echo "Running second BaseRecalibrator pass for ${id}"
    echo ""

    gatk BaseRecalibrator -R "${params.reference}" -I "${recal_bam}" --known-sites "${params.KnownSites}" -O "${id}_${params.Ref_Abbr}_Recal_Data.after.table"

    # Generate the plot

    gatk AnalyzeCovariates -before "${data_table}" -after "${id}_${params.Ref_Abbr}_Recal_Data.after.table" -csv "${id}_Recalibration.csv"

    echo ""
    echo "Done with ${id}"
    echo ""

    """
    }

process HapCall {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.VCF}", mode: 'copy'

    input:
    tuple val(id), path(recal_marked_bam), val(chrom)
    path params.reference

    output: 
    tuple val(id), val(chrom),
          path("${id}.${chrom}.${params.Ref_Abbr}.g.vcf.gz"),
          path("${id}.${chrom}.${params.Ref_Abbr}.g.vcf.gz.tbi")

    script:
    """
    echo "Job info:"
    echo "- Processing sample: ${id}"
    echo "- Processing chromosome: ${chrom}"

    echo "\$(date) - Starting HaplotypeCaller for ${id} on chromosome ${chrom}"

    gatk --java-options "-Xmx100g -XX:+UseParallelGC" HaplotypeCaller --native-pair-hmm-threads 64 -R "${params.reference}" -I "${recal_marked_bam}" -L "${chrom}" -O "${id}.${chrom}.${params.Ref_Abbr}.g.vcf.gz" -ERC GVCF

    echo "\$(date) - Finished HaplotypeCaller on ${id} chromosome ${chrom}"
    """
    }

process chr_GenomicDB {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.DB}", mode: 'copy'

    input:
    tuple val(chrom), path(gVCFs)

    output: 
    tuple val(chrom), path("${chrom}_DB")

    script:
    """
    echo "\$(date) - Starting GenomicsDBImport for chromosome ${chrom}"

    CHR_DB_DIR="${chrom}_DB"

    GVCF_ARGS=""
    for GVCF_PATH in ${gVCFs}; do
        if [[ "\$GVCF_PATH" == *.g.vcf.gz ]]; then
            GVCF_ARGS="\$GVCF_ARGS -V \$GVCF_PATH"
        fi
    done

    # Run GenomicsDBImport for this chromosome only
    gatk --java-options "-Xmx100g -XX:+UseParallelGC" GenomicsDBImport --reader-threads 128 \$GVCF_ARGS --genomicsdb-workspace-path "\$CHR_DB_DIR" -L "${chrom}"

    echo "\$(date) - Finished GenomicsDBImport for chromosome ${chrom}"
    """
    }

process genotype {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.final_VCF}", mode: 'copy'

    input:
    tuple val(chrom), path(DB)
    path params.reference

    output: 
    tuple val(chrom), path("${params.vcf_prefix}.${chrom}.${params.Ref_Abbr}.raw.vcf.gz"),
          path("${params.vcf_prefix}.${chrom}.${params.Ref_Abbr}.raw.vcf.gz.tbi")

    script:
    """
    echo "\$(date) - Starting Genotyping GVCFs for ${chrom}"

    gatk --java-options "-Xmx80g -XX:+UseParallelGC" GenotypeGVCFs -R "${params.reference}" -V "gendb://${DB}" -O "${params.vcf_prefix}.${chrom}.${params.Ref_Abbr}.raw.vcf.gz"

    echo "\$(date) - Finished Genotyping GCVFs for ${chrom}"
    """
    }

process MergeVcfs {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.final_VCF}", mode: 'copy'

    input:
    path(chr_vcfs)

    output: 
    tuple path("${params.vcf_prefix}_whole_genome.${params.Ref_Abbr}.raw.vcf.gz"),
          path("${params.vcf_prefix}_whole_genome.${params.Ref_Abbr}.raw.vcf.gz.tbi")


    script:
    """
    echo "\$(date) - Starting Collating VCFs"

    VCF_ARGS=""
        for vcf in ${chr_vcfs}; do
            VCF_ARGS="\$VCF_ARGS -I \$vcf"
        done


    gatk MergeVcfs \$VCF_ARGS -O "${params.vcf_prefix}_whole_genome.${params.Ref_Abbr}.raw.vcf.gz"

    echo "Done!"
    """
    }

workflow {
    reads = Channel.fromFilePairs("${params.reads_dir}/*_{${params.read1_tag},${params.read2_tag}}.fq.gz", flat: true)
    clean_reads = fastp(reads).map { id, r1, r2, qc -> tuple(id, r1, r2) }
    ref_index = file(params.Ref_index + ".*")
    mapped = mapping(clean_reads, ref_index)
    marked_bam = DupValIndx(mapped).map { id, marked, index, metrics, validation ->  tuple(id, marked)}
    reference_ch = Channel.value(file("${baseDir}/Genome/Reference"))
    known_sites_ch = Channel.value(file("${baseDir}/Genome"))
    marked_bam
          .combine(reference_ch)
          .combine(known_sites_ch)
          .set { recal_inputs }
    recal_bam = recalibrate(recal_inputs).map { id, recal, index, data, val -> tuple(id, recal, data) }
    recal_bam
          .combine(reference_ch)
          .combine(known_sites_ch)
          .set { BQSR_inputs }
    BQSR(BQSR_inputs)
    chromosomes = Channel
          .fromPath(params.chr_file)
          .splitText()
          .map { it.trim() }
    recal_marked_bam = recal_bam.map { id, recal, data -> tuple(id, recal) }
    hapcall_inputs = recal_marked_bam.combine(chromosomes)
    DB_inputs = HapCall(hapcall_inputs, file(params.reference))
          .map { id, chrom, gvcf, tbi -> tuple(chrom, [gvcf, tbi]) }  
          .groupTuple(by: 0)
          .map { chrom, files_list -> tuple(chrom, files_list.flatten()) }
    database = chr_GenomicDB(DB_inputs)
    chr_vcfs = genotype(database, file(params.reference)).map { chrom, vcf, index -> vcf }
    chr_vcfs.collect().set { all_vcfs }
    MergeVcfs(all_vcfs)
}
