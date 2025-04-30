#!/bin/bash
#SBATCH --job-name=agat-ahya-convert
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o agat-convert-ahya-%j.out
#SBATCH -e agat-convert-ahya-%j.err

module load uri/main
module load conda/latest
source $(conda info --base)/etc/profile.d/conda.sh
conda activate agat-env
conda install -c bioconda agat

# Optional: check if AGAT is really there
which agat_convert_sp_gff2gbk.pl

agat_convert_sp_gff2gbk.pl \
  --gff Ahyacinthus.coding.gff3 \
  --fasta Ahyacinthus.proteins.fasta \
  --output ahya.gb

echo "AGAT conversion complete: ahya.gb created."
