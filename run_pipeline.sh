#!/bin/bash

# SAILOR Nextflow Pipeline Execution Script
# Modify the parameters below according to your setup

#=============================================================================
# JOB SCHEDULER CONFIGURATION (SLURM example)
# Uncomment and modify if using SLURM
#=============================================================================
#SBATCH --job-name=SAILOR_Pipeline
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=./logs/sailor_pipeline_%j.stdout
#SBATCH --error=./logs/sailor_pipeline_%j.stderr

# Activate conda environment
# source ~/.bashrc
# conda activate sailor-pipeline

#=============================================================================
# PIPELINE PARAMETERS - MODIFY THESE PATHS
#=============================================================================

# Genome and reference files
GENOME_DIR="/path/to/your/genome_directory"
GENOME="${GENOME_DIR}/genome.fa"
GENOME_INDEX="${GENOME_DIR}/genome.fa.fai"
GENOME_DIR_STAR="${GENOME_DIR}/star_index"
GTF="${GENOME_DIR}/annotation.gtf"
KNOWN_SNP="${GENOME_DIR}/known_snps.bed"

# Data directories
INPUT_DIR="/path/to/your/fastq_files"
OUTPUT_DIR="/path/to/your/output_directory"
WORK_DIR="${OUTPUT_DIR}/nextflow_work"

# Sample sheets
SAMPLE_SHEET="assets/sample_sheet.csv"
DESEQ2_SHEET="assets/sample_sheet_deseq2.csv"

# Analysis parameters
GENOME_BUILD="hg38"
EDITOR="APOBEC"
PROJECT_NAME="My_SAILOR_Project"
PADJ_THRESHOLD="0.05"
DIFF_THRESHOLD_DESEQ2="1"

# Workflow flags
WORKFLOW_DE="DESEQ2"        # Set to "NULL" to skip
WORKFLOW_SAILOR="SAILOR"    # Set to "NULL" to skip

# Execution profile
PROFILE="standard"          # Options: standard, slurm, docker, singularity

#=============================================================================
# FILE EXISTENCE CHECKS
#=============================================================================
echo "Checking required files..."

check_file() {
    if [ ! -f "$1" ] && [ ! -d "$1" ]; then
        echo "ERROR: Required file/directory does not exist: $1"
        exit 1
    else
        echo "✓ Found: $1"
    fi
}

check_file "$GENOME"
check_file "$GENOME_INDEX"
check_file "$GENOME_DIR_STAR"
check_file "$GTF"
check_file "$KNOWN_SNP"
check_file "$SAMPLE_SHEET"
check_file "$DESEQ2_SHEET"

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

echo "All required files found!"

#=============================================================================
# CREATE OUTPUT DIRECTORIES
#=============================================================================
mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$WORK_DIR")"
mkdir -p logs

#=============================================================================
# RUN PIPELINE
#=============================================================================
echo "Starting SAILOR Nextflow Pipeline..."
echo "Output directory: $OUTPUT_DIR"
echo "Work directory: $WORK_DIR"

nextflow run main.nf \
    --genome_dir "$GENOME_DIR" \
    --genome "$GENOME" \
    --genome_index "$GENOME_INDEX" \
    --genomeDir "$GENOME_DIR_STAR" \
    --gtf "$GTF" \
    --known_snp "$KNOWN_SNP" \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --sample_sheet "$SAMPLE_SHEET" \
    --deseq2_test_sheet "$DESEQ2_SHEET" \
    --genome_build "$GENOME_BUILD" \
    --editor "$EDITOR" \
    --project_name "$PROJECT_NAME" \
    --padj_threshold "$PADJ_THRESHOLD" \
    --diff_threshold_deseq2 "$DIFF_THRESHOLD_DESEQ2" \
    --workflow_DE "$WORKFLOW_DE" \
    --workflow_SAILOR "$WORKFLOW_SAILOR" \
    -profile "$PROFILE" \
    -work-dir "$WORK_DIR" \
    -resume

echo "Pipeline execution completed!"
echo "Check the output directory for results: $OUTPUT_DIR"
echo "Execution reports available at:"
echo "  - Timeline: ${OUTPUT_DIR}/timeline.html"
echo "  - Report: ${OUTPUT_DIR}/report.html"
echo "  - Trace: ${OUTPUT_DIR}/trace.txt"