#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=ATAC-seq_alignment_and_file_conversion
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job. Replace <PI GROUP> with your own group.
#SBATCH --account=your_account_name
### REQUIRED. Set the partition for your job. This is a job queue
#SBATCH --partition=standard
### REQUIRED. Set the number of nodes
#SBATCH --nodes=1
### REQUIRED. Set the number of CPUs that will be used for this job. 
#SBATCH --ntasks=12
### REQUIRED. Set the memory required for this job.
#SBATCH --mem-per-cpu=5gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=80:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

# Activate conda environment
source ~/.bashrc
eval "$(micromamba shell hook --shell=bash)"
micromamba activate ATAC-seq_env

# Load required modules
# For example
module load bowtie2/2.5.4
module load samtools/1.19.2

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
run_if_exists ./path/to/your/bash_file/atac-seq_mouse.sh

echo "========== Job finished at: $(date) =========="

