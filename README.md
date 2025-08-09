# SAILOR Nextflow Pipeline

A Nextflow pipeline for RNA editing analysis using SAILOR (Systematic Analysis of Inosine, Logging Operations, and Reviewing) methodology. This pipeline integrates RNA-seq mapping, RNA editing detection, and differential expression analysis.

## Overview

This pipeline performs:
- **Pre-processing**: FastQ quality control and trimming
- **Alignment**: RNA-seq reads mapping with STAR aligner
- **Quality Control**: BAM file quality assessment with Qualimap
- **RNA Editing Analysis**: SAILOR-based editing site detection
- **Differential Expression**: DESeq2-based differential gene expression analysis
- **Visualization**: Generation of plots and reports

## Requirements

### Software Dependencies
- [Nextflow](https://nextflow.io) (>= 20.04.0)
- [Conda/Mamba](https://conda.io) or [Docker](https://docker.com)/[Singularity](https://sylabs.io/singularity)

### Reference Data Required
- Reference genome FASTA file
- Reference genome index (.fai)
- STAR genome index directory
- Gene annotation GTF file
- Known SNPs BED file (for filtering)

## Installation

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd sailor-nextflow-pipeline
   ```

2. **Install dependencies using conda/mamba:**
   
   **Option A: Full environment (recommended)**
   ```bash
   # Create environment from YAML file
   conda env create -f environment.yml
   conda activate sailor-pipeline
   ```
   
   **Option B: Minimal environment (faster install)**
   ```bash
   # Create minimal environment
   conda env create -f environment-minimal.yml
   conda activate sailor-pipeline-minimal
   ```
   
   **Option C: Manual installation**
   ```bash
   # Create and activate environment
   conda create -n sailor-pipeline python=3.9
   conda activate sailor-pipeline
   
   # Install core tools
   conda install -c bioconda nextflow star samtools fastp qualimap multiqc htseq
   conda install -c bioconda gatk4 picard vcftools bedtools
   conda install -c conda-forge r-base r-argparse
   conda install -c bioconda bioconductor-deseq2 r-ggplot2 r-dplyr
   ```

3. **Or use container systems:**
   - Docker: Use `-profile docker`
   - Singularity: Use `-profile singularity`

4. **Verify installation:**
   ```bash
   # Test Nextflow
   nextflow -version
   
   # Test key tools
   STAR --version
   samtools --version
   fastp --version
   Rscript -e "library(DESeq2); sessionInfo()"
   
   # Test pipeline syntax
   nextflow run main.nf --help
   ```

## Usage

### 1. Prepare Input Files

Create sample sheets in the `assets/` directory:

**Sample Sheet (`assets/sample_sheet.csv`):**
```csv
Group,Sample,File_name,Folder_name,fastq_1,fastq_2
Control,Control_1,C1,Control_samples,C1_R1.fastq.gz,C1_R2.fastq.gz
Control,Control_2,C2,Control_samples,C2_R1.fastq.gz,C2_R2.fastq.gz
Treatment,Treatment_1,T1,Treatment_samples,T1_R1.fastq.gz,T1_R2.fastq.gz
Treatment,Treatment_2,T2,Treatment_samples,T2_R1.fastq.gz,T2_R2.fastq.gz
```

**DESeq2 Comparison Sheet (`assets/sample_sheet_deseq2.csv`):**
```csv
Condition,Control
Treatment,Control
```

### 2. Run the Pipeline

**Basic usage:**
```bash
nextflow run main.nf \
  --genome_dir /path/to/genome \
  --genome /path/to/genome/genome.fa \
  --genome_index /path/to/genome/genome.fa.fai \
  --genomeDir /path/to/star_index \
  --gtf /path/to/annotation.gtf \
  --known_snp /path/to/known_snps.bed \
  --input_dir /path/to/fastq_files \
  --output_dir /path/to/results \
  --sample_sheet assets/sample_sheet.csv \
  --deseq2_test_sheet assets/sample_sheet_deseq2.csv
```

**With SLURM profile:**
```bash
nextflow run main.nf \
  --genome_dir /path/to/genome \
  --genome /path/to/genome/genome.fa \
  --genome_index /path/to/genome/genome.fa.fai \
  --genomeDir /path/to/star_index \
  --gtf /path/to/annotation.gtf \
  --known_snp /path/to/known_snps.bed \
  --input_dir /path/to/fastq_files \
  --output_dir /path/to/results \
  --sample_sheet assets/sample_sheet.csv \
  --deseq2_test_sheet assets/sample_sheet_deseq2.csv \
  -profile slurm
```

### 3. Pipeline Parameters

#### Required Parameters:
- `--genome_dir`: Directory containing reference genome files
- `--genome`: Reference genome FASTA file
- `--genome_index`: Reference genome index (.fai file)  
- `--genomeDir`: STAR genome index directory
- `--gtf`: Gene annotation GTF file
- `--known_snp`: Known SNPs BED file
- `--input_dir`: Directory containing input FASTQ files
- `--sample_sheet`: Sample information CSV file

#### Optional Parameters:
- `--output_dir`: Output directory (default: './results')
- `--genome_build`: Genome build version (default: 'hg38')
- `--editor`: RNA editor type (default: 'APOBEC')
- `--workflow_DE`: Run DESeq2 analysis ('DESEQ2' or 'NULL', default: 'DESEQ2')
- `--workflow_SAILOR`: Run SAILOR analysis ('SAILOR' or 'NULL', default: 'SAILOR')
- `--diff_threshold_deseq2`: Log2 fold change threshold for DESeq2 (default: 1)
- `--padj_threshold`: Adjusted p-value threshold (default: 0.05)
- `--project_name`: Project name for outputs (default: 'SAILOR_analysis')

### 4. Execution Profiles

- `standard`: Local execution (default)
- `slurm`: SLURM cluster execution
- `docker`: Docker container execution
- `singularity`: Singularity container execution

## Output Structure

```
results/
├── timeline.html          # Nextflow execution timeline
├── report.html           # Nextflow execution report
├── trace.txt             # Nextflow trace file
├── dag.svg               # Pipeline DAG
├── bam_files/            # Merged BAM files
├── sailor/               # SAILOR editing analysis results
├── deseq2/               # DESeq2 differential expression results
└── plots/                # Generated visualizations
```

## Workflow Details

### 1. Pre-processing (CLEAN_FASTQ)
- Quality trimming with fastp
- Adapter removal
- Quality control metrics

### 2. RNA-seq Mapping (RNASEQ_MAPPING_STAR)
- 2-pass STAR alignment
- Unique read filtering
- BAM indexing

### 3. Quality Assessment (RUN_QUALIMAP)
- BAM file quality metrics
- Coverage analysis
- Insert size distribution

### 4. BAM Merging (MERGE_BAMS)
- Merge BAM files by sample
- Re-indexing merged files

### 5. RNA Editing Analysis (SAILOR)
- Site-specific editing detection
- Strand-aware analysis
- Coverage filtering

### 6. Differential Expression (DESEQ2)
- Gene count quantification with HTSeq
- Differential expression analysis
- Visualization (PCA, volcano plots)

## Troubleshooting

### Common Issues:

1. **Memory errors**: Increase memory allocation in `nextflow.config`
2. **File not found**: Check all file paths are absolute and accessible
3. **Permission denied**: Ensure write permissions for output directory
4. **STAR index errors**: Verify STAR index was built with compatible version

### Getting Help:

1. Check Nextflow execution reports (`timeline.html`, `report.html`)
2. Review log files in the work directory
3. Verify input file formats match expected schemas

## Citation

If you use this pipeline, please cite:
- **SAILOR**: Washbourne et al. (original SAILOR methodology)
- **Nextflow**: Di Tommaso et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017)

## License

This project is licensed under the Mozilla Public License 2.0 - see the original license headers in the source files for details.