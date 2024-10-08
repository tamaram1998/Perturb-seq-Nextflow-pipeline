# Perturb-Seq Pipeline

Perturb-Seq Pipeline is designed for processing Perturb-Seq data, encompassing both primary and secondary analysis stages. It automates the workflow of quality control, read alignment, and gene quantification, providing a streamlined solution for genomic data processing.

For additional information on proposed architecture on AWS, refer to [this documentation](aws_architecture/architecture.md).

## Pipeline Components

This pipeline is designed to handle the primary processing steps for Perturb-seq data, a high-throughput method that combines CRISPR-based perturbations with single-cell RNA sequencing (scRNA-seq).

Perturb-seq, also known as CROP-seq or CRISP-seq, allows for the simultaneous measurement of gene expression and CRISPR-induced genetic perturbations at the single-cell level. This powerful technique enables researchers to study gene function and regulatory networks at scale.

While our pipeline doesn't include specialized steps for guide RNA (gRNA) quantification or cell barcode processing, it provides essential preprocessing for Perturb-seq data.

This Perturb-Seq Pipeline handles the initial processing of Perturb-seq data, preparing it for subsequent specialized analyses. Note that additional steps are required for:
- Identification and quantification of guide RNAs
- Cell barcode and UMI (Unique Molecular Identifier) processing
- Integration of perturbation information with gene expression data

## Prerequisites

- Docker
- Nextflow (version 20.04.0 or later)
- Java (version 8 or later)

## Installation

This pipeline is containerized using Docker to ensure consistency and reproducibility across different environments.

Clone this repository:

    git clone https://github.com/yourusername/genexomics-pipeline.git

    cd genexomics-pipeline

Build the Docker image:

    docker build -t genexomics-pipeline:latest .

Running the Docker image:

    docker run -it genexomics-pipeline:latest

## Usage

Ensure your input data (FASTQ files, reference genome and GTF file) are in a known location on your system.

Example pipeline command with custom parameters:

- `--reads`: Path to input reads (FASTQ format). Use wildcard patterns to specify paired-end reads
- `--genome`: Path to reference genome (FASTA format)
- `--gtf`: Path to gene annotation file (GTF format)
- `--outdir`: Output directory

    ```
    nextflow run main.nf genexomics-pipeline:latest \
        --reads "test_data/reads_test.fastq.gz" \
        --genome "test_data/genome_test.fa" \
        --gtf "test_data/gtf_test.gtf" \
        --outdir "genexomics-results"
    ```

## Pipeline Stages

### 1. Quality Control (FASTQC)

FastQC performs quality control checks on raw sequence data. It provides an overview of potential quality issues in your data.

**Input**: Raw FASTQ files
**Output**: HTML reports and ZIP archives containing QC results

### 2. Genome Indexing (STAR_INDEX)

This step creates an index of the reference genome, which is required for the alignment process.

**Input**: Reference genome FASTA file and GTF annotation
**Output**: STAR index files

### 3. Read Alignment (ALIGN)

STAR aligner maps the sequencing reads to the reference genome.

**Input**: FASTQ files and STAR index
**Output**: BAM files containing aligned reads

### 4. Gene Quantification (FEATURECOUNTS)

featureCounts quantifies mapped reads for genomic features such as genes, exons, promoters, etc.

**Input**: BAM files and GTF annotation
**Output**: Tab-delimited text files with read counts per feature

### 5. Report Generation (MULTIQC)

MultiQC aggregates results from the previous steps into a single HTML report.

**Input**: Output from FastQC, STAR, and featureCounts
**Output**: HTML report summarizing the analysis results

## Output

The pipeline produces the following outputs:

1. FastQC

    FastQC generates HTML reports for each input FASTQ file. These reports provide insights into the quality of your sequencing data.

    [Example output](example_outputs/fastqc_report.html)

2. STAR Alignment 

    STAR generates log files that provide information about the alignment process.

    [Example output](example_outputs/alignment.sortedByCoord.out.bam)

3. featureCounts 

    featureCounts generates a tab-delimited file containing gene-level quantification data.

    [Example output](example_outputs/featurecounts_result.txt)

4. MultiQC Report

    MultiQC generates a comprehensive HTML report that aggregates results from all other tools used in the pipeline.

    [Example output](example_outputs/multiqc_report.html)


## Test Dataset

The test samples required for this pipeline are stored in an Amazon S3 bucket due to their size. There are three main input files:

1. Genome test file (FASTA format):
   - Filename: genome_test.fa
   - Size: 237.7 MB
   - [Download here](https://pipeline-test-datasets.s3.amazonaws.com/test_datasets/genome_test.fa)

2. GTF test file:
   - Filename: gtf_test.gtf
   - Size: 40.5 MB
   - [Download here](https://pipeline-test-datasets.s3.amazonaws.com/test_datasets/gtf_test.gtf)

3. Small FASTQ test file:
   - Filename: SRR23584814_small.fastq
   - Size: 235.1 KB
   - [Download here](https://pipeline-test-datasets.s3.amazonaws.com/test_datasets/SRR23584814_small.fastq)


To download all three files together, run the following command:

```bash
./download_data.sh
```
This command will automatically retrieve all the necessary test files in one go. Make sure the script is executable and your environment is properly configured to download the files.


## References

FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

STAR: https://github.com/alexdobin/STAR

featureCounts: http://subread.sourceforge.net/

MultiQC: https://multiqc.info/

Nextflow: https://www.nextflow.io/
