#!/bin/bash
#SBATCH --job-name=toxcodan-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o toxcodan-pever-%j.out
#SBATCH -e toxcodan-pever-%j.error

module load conda
conda activate ToxcodanGenome

