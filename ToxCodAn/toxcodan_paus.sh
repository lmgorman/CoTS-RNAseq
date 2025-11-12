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

# Ensure Conda uses /work directories
export CONDA_PKGS_DIRS=/work/pi_hputnam_uri_edu/lgorman/.conda/pkgs
export CONDA_ENVS_DIRS=/work/pi_hputnam_uri_edu/lgorman/.conda/envs

# -----------------------------
# Load modules
# -----------------------------
module load uri/main
module load conda/latest

# Activate environment
conda activate ToxcodanGenome


# Run script
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py \
    -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna \
    -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta \
    -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results \
    -c 10
