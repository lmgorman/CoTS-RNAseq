#!/bin/bash
#SBATCH --job-name=STAR_genome_index_Ahya
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1

# Load modules (add gffread module if needed)
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0
module load gffread/0.12.7

echo "[$(date)] Starting STAR genome index build with GFF3 annotation"

# Variables
GENOME_FA="GCA_964291705.1_jaAcrHyac4.1_genomic.fna"
GTF_FILE="GCA_964291705.1_jaAcrHyac4.1_genomic.gff3"
INDEX_DIR="star_index_jaAcrHyac4.1"

# Create output directory
mkdir -p "$INDEX_DIR"

# Run STAR genomeGenerate with annotation
STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF_FILE" \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 99

echo "[$(date)] Genome index build completed"
