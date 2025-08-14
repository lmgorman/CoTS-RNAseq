#!/bin/bash
#SBATCH --job-name=STAR_genome_index_Ahya
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH -p uri-cpu
#SBATCH --time=24:00:00
#SBATCH -o slurm-%j.out
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1

# Load modules
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0
module load gffread/0.12.7

# Set custom scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/star_index_jaAcrHyac4.1_$SLURM_JOB_ID
mkdir -p "$SCRATCHDIR"

echo "[$(date)] Converting GFF3 to GTF"
gffread /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff.gff \
  -T -o GCA_964291705.1_jaAcrHyac4.1_genomic.gtf

echo "[$(date)] Starting STAR genome index build in scratch: $SCRATCHDIR"

# Variables
GENOME_FA="GCA_964291705.1_jaAcrHyac4.1_genomic.fna"
GTF_FILE="GCA_964291705.1_jaAcrHyac4.1_genomic.gtf"
INDEX_DIR="$SCRATCHDIR"

# Run STAR genomeGenerate
STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF_FILE" \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 99

# Copy STAR index back to permanent storage
FINAL_OUTPUT_DIR="/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/star_index_jaAcrHyac4.1"
mkdir -p "$FINAL_OUTPUT_DIR"
cp -r "$INDEX_DIR"/* "$FINAL_OUTPUT_DIR"

echo "[$(date)] Genome index build completed and copied to $FINAL_OUTPUT_DIR"

