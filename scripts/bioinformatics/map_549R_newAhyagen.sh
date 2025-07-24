#!/bin/bash
#SBATCH --job-name=Ahya_align_star
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH --time=24:00:00
#SBATCH -o slurm-%j.out
#SBATCH -D /work/pi_hputnam_uri_edu/20250107_COTS_LG/

# Load modules
echo "Loading programs" $(date)
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0

# Input/output
R1="/work/pi_hputnam_uri_edu/20250107_COTS_LG/fastp_trimmed/Ahya/549R_R1_001_trimmed.fastq.gz"
R2="/work/pi_hputnam_uri_edu/20250107_COTS_LG/fastp_trimmed/Ahya/549R_R2_001_trimmed.fastq.gz"
GENOME_DIR="/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1"
OUT_PREFIX="Ahya_STAR_new"

# Run STAR
echo "Starting read alignment." $(date)

STAR --runMode alignReads \
--genomeDir "$GENOME_DIR" \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn "$R1" "$R2" \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix "$OUT_PREFIX"

echo "Alignment of Trimmed Seq data complete." $(date)
