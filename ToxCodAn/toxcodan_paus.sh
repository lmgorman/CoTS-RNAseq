#!/bin/bash
#SBATCH --job-name=toxcodan-paus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/logs/toxcodan-paus-%j.out
#SBATCH -e /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/logs/toxcodan-paus-%j.error

# -----------------------------
# Load modules
# -----------------------------
module load uri/main
module load conda/latest

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home/lucy_gorman_uri_edu/.conda/envs/ToxcodanGenome


# Run script
/home/lucy_gorman_uri_edu/.conda/envs/ToxcodanGenome/bin/python \
    /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py \ 
    -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna \
    -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta.fasta \
    -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results \
    -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta \
    -c 10
