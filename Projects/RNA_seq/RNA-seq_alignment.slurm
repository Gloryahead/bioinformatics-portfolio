#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources to run your job
# --------------------------------------------------------------

### Optional. Set the job name
#SBATCH --job-name=RNA_seq_STAR_align
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
#SBATCH --time=22:00:00

# --------------------------------------------------------------
### PART 2: Setup and run your job
# --------------------------------------------------------------

echo "========== Job started at: $(date) =========="

# Activate conda environment
source ~/.bashrc
eval "$(micromamba shell hook --shell=bash)"
micromamba activate RNA-seq

# Load required modules
module load star/2.7.11b
module load samtools/1.22.1
module load subread/2.1.1

# Run bash script
#bash RNA_seq/RNA-seq_star_align_samtools.sh
bash RNA_seq/RNA-seq_count_transcript.sh

echo "========== Job finished at: $(date) =========="

