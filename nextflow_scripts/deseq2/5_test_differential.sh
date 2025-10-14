#!/bin/bash
#SBATCH --job-name=5_test_differential
#SBATCH --time=8:00:00
#SBATCH --array=1-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=./logs/5_test_differential_%j.stdout
#SBATCH --error=./logs/5_test_differential_%j.stderr
# initialize conda environment

source ~/.bashrc
conda activate nextflow_mamba

# comparisons=(
#     "shCon_Ribo,shCon_MIB"
#     "shcon_Ribo,shcon_MIB"
#     "shMOV10_Ribo,shMOV10_MIB"
#     "shMSI2_Ribo,shMSI2_MIB"
#     "shSYN_Ribo,shSYN_MIB"
#     "shMSI2_Ribo,shcon_Ribo"
#     "shMOV10_Ribo,shCon_Ribo"
#     "shSYN_Ribo,shCon_Ribo"
#     "shMSI2_MIB,shcon_MIB"
#     "shMOV10_MIB,shCon_MIB"
#     "shSYN_MIB,shCon_MIB"
#     "shCon_Ribo,shcon_Ribo"
#     "shMSI2_Ribo,shMOV10_Ribo"
#     "shMSI2_Ribo,shSYN_Ribo"
#     "shSYN_Ribo,shMOV10_Ribo"
#     "shMSI2_MIB,shCon_MIB"
#     "shCon_MIB,shcon_MIB"
# )


comparisons=(
    "shCon_Ribo,shCon_MIB"
    "shcon_MIB,shCON_MIG"
    "shCon_MIB,shCON_MIG"
    "shcon_MIB,shCON_MSI2_ADAR"
    "shCon_MIB,shCON_MSI2_ADAR"    
)


# Get the SLURM array task ID
index=${SLURM_ARRAY_TASK_ID}

# Parse ctrl and cond from the comparisons array
IFS=',' read -r cond ctrl  <<< "${comparisons[$index]}"

# R: format output
# Rscript 5_test_differential.R --counts_file /data1/kharasm/baalii/RiboSTAMP_MOLM/output_data/sailor_new/deseq2/4_formatted_output_counts.csv  --ctrl $ctrl --cond $cond

echo "Running DESeq2 for comparison: $cond vs $ctrl"

Rscript 5_test_differential.R --counts_file /data1/kharasm/baalii/RiboSTAMP_MOLM/output_data/sailor_new/deseq2/4_formatted_output_ribostamp_hypertribe_controls_counts.csv --ctrl $ctrl --cond $cond


