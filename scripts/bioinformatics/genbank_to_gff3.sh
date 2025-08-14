#!/bin/bash
#SBATCH --job-name=genbank_to_gff3_Ahya
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman

# Clone github repo
git clone https://github.com/bioperl/bioperl-live.git

#Run bioperl
module load uri/main
module load BioPerl/1.7.8-GCCcore-13.3.0
perl /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/bin/bp_genbank2gff3 /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff
