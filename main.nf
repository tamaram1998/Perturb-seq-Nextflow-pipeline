#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Input parameters
params.reads = "/data/other_data/SRR23584814_small_fixed.fastq"
params.genome = "/data/reference_data/genome_test.fa"
params.gtf = "/data/reference_data/gtf_test.gtf"
params.outdir = "/data/genexomics-results"

// Validate input parameters
if (!params.reads || !params.genome || !params.gtf) {
    error "Missing required parameters. Please provide --reads, --genome, and --gtf."
}

// Create input channels with validation
reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
                  .map { file -> tuple(file.simpleName, file) }
genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)

// Log parameter values
log.info """
    ===================================
    genexomics-pipeline
    ===================================
    reads        : ${params.reads}
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    outdir       : ${params.outdir}
    """
    .stripIndent()


// Processes
process VALIDATE_FASTQ {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("validated_${reads}")

    script:
    """
    zcat $reads | head -n 400 | awk 'NR%4==0{printf "%s\\n", \$0}' | grep -P '^[!-~]+\$' > /dev/null && cp $reads validated_${reads} || echo "Invalid FASTQ file" >&2
    """
}

process STAR_INDEX {
    tag "STAR_INDEX"
    publishDir "${params.outdir}/star_index", mode: 'copy'
    
    input:
    path genome
    path gtf

    output:
    path "star_index", type: 'dir'

    script:
    """
    mkdir star_index
    
    if [[ $genome == *.gz ]]; then
        gunzip -c $genome > genome.fa
        GENOME_FILE=genome.fa
    else
        GENOME_FILE=$genome
    fi

    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles \$GENOME_FILE \
         --outFileNamePrefix star_index/ \
         --sjdbGTFfile $gtf \
         --genomeSAindexNbases 12 \
         --runThreadN 4

    if [[ $genome == *.gz ]]; then
        rm genome.fa
    fi
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
 
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}_fastqc_out"

    script:
    """
    mkdir ${sample_id}_fastqc_out
    fastqc -o ${sample_id}_fastqc_out -f fastq -q ${reads} || echo "FastQC failed for ${reads}" >&2
    """
}

process ALIGN {
    tag "STAR on $sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path star_index
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR --runThreadN 4 \
         --genomeDir $star_index \
         --readFilesIn ${reads} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_ \
         --outSAMtype BAM SortedByCoordinate \
         --genomeSAindexNbases 12 \
         --sjdbGTFfile ${gtf}
    """
}

process FEATURECOUNTS {
    tag "featureCounts on $sample_id"
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path gtf
    
    output:
    path "${sample_id}_counts.txt"

    script:
    """
    featureCounts -a ${gtf} -o ${sample_id}_counts.txt ${bam}
    """
}

process MULTIQC {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path '*'
    path 'aligned/*'
    path 'counts/*'
    
    
    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

// Workflow
workflow {
    fastqc_ch = FASTQC(reads_ch)
    star_index_ch = STAR_INDEX(genome_ch, gtf_ch)
    aligned_ch = ALIGN(reads_ch, star_index_ch, gtf_ch)
    counts_ch = FEATURECOUNTS(aligned_ch, gtf_ch)
    
    MULTIQC(
        fastqc_ch.collect(),
        aligned_ch.map { it[1] }.collect(),
        counts_ch.collect()
    )
}