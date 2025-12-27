#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=CR_alignment_and_file_conversion
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job. Replace <PI GROUP> with your own group.
#SBATCH --account=clsmith1
### REQUIRED. Set the partition for your job. This is a job queue
#SBATCH --partition=standard
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the number of CPUs that will be used for this job. 
#SBATCH --ntasks=8
### REQUIRED. Set the memory required for this job.
#SBATCH --mem-per-cpu=5gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=160:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

# Activate conda environment
source ~/.bashrc
eval "$(micromamba shell hook --shell=bash)"
micromamba activate Your_env

# Load required modules
module load bowtie2/2.5.4
module load samtools/1.19.2
module load bedtools/2.31.1
module load deeptools/3.5.4
module load openjdk/22.0.1
module load macs2/2.2.9.1
module load homer/5.1

# Make all scripts executable
chmod u+x /path/to/your/bash_file/*.sh

# Function to run script if it exists
run_if_exists() {
    SCRIPT_PATH=$1
    if [[ -f "$SCRIPT_PATH" ]]; then
        echo "Running: $SCRIPT_PATH"
        bash "$SCRIPT_PATH"
    else
        echo "WARNING: Script not found: $SCRIPT_PATH"
    fi
}

# Run processing scripts
run_if_exists ./path/to/your/bash_file/run_bowtie_file_conversion_mark_remove_dup.sh
run_if_exists ./path/to/your/bash_file/bigwig.sh
run_if_exists ./path/to/your/bash_file/bigwig_RPKM.sh
run_if_exists ./path/to/your/bash_file/homer.sh
run_if_exists ./path/to/your/bash_file/run_macs2_peaks.sh
run_if_exists ./path/to/your/bash_file/heatmaps.sh

echo "========== Job finished at: $(date) =========="

