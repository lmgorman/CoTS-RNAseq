#!/bin/bash
#SBATCH --job-name=Ahya_align_star
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH -D work/pi_hputnam_uri_edu/20250107_COTS_LG/

#load modules
echo "Loading programs" $(date)
module load uri/main
module load module load genometools/1.6.5

