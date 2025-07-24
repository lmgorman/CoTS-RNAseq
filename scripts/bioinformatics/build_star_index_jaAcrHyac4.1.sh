#!/bin/bash
#SBATCH --job-name=build_star_index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH -p short
#SBATCH -o slurm-%j-build-index.out
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/

# Load modules (add gffread module if needed)
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0
module load gffread/0.12.7

echo "[$(date)] Starting STAR genome index build with GFF3 annotation"

# Variables
GENOME_FA_GZ="GCA_964291705.1_jaAcrHyac4.1_genomic.fna.gz"
GENOME_FA="genome.fa"
GFF3_FILE="GCA_964291705.1_jaAcrHyac4.1_genomic.gff3"   
GTF_FILE="annotation.gtf"
INDEX_DIR="star_index_jaAcrHyac4.1"

# Uncompress genome FASTA if needed
if [ ! -f "$GENOME_FA" ]; then
    echo "Uncompressing genome FASTA..."
    gunzip -c "$GENOME_FA_GZ" > "$GENOME_FA"
fi

# Convert GFF3 to GTF
if [ ! -f "$GTF_FILE" ]; then
    echo "Converting GFF3 to GTF..."
    gffread "$GFF3_FILE" -T -o "$GTF_FILE"
    if [ $? -ne 0 ]; then
        echo "ERROR: gffread failed to convert GFF3 to GTF"
        exit 1
    fi
else
    echo "GTF file already exists, skipping conversion."
fi

# Create output directory
mkdir -p "$INDEX_DIR"

# Run STAR genomeGenerate with annotation
STAR --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF_FILE" \
    --runThreadN 8 \
    --sjdbOverhang 99

echo "[$(date)] Genome index build completed"
