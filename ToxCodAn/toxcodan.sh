#!/bin/bash
#SBATCH --job-name=toxcodan-paus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o toxcodan-paus-%j.out
#SBATCH -e toxcodan-paus-%j.error

#Load modules
module load uri/main
module load conda/latest
conda activate ToxcodanGenome

#Run script
python toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta

