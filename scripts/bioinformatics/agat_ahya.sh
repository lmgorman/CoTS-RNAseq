#!/bin/bash
#SBATCH --job-name=agat-ahya
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o agat-ahya-%j.out
#SBATCH -e agat-ahya-%j.error

set -e  # Exit on first error

# Define scratch directory
SCRATCHDIR=/scratch/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR

echo "[$(date)] Job started in $SCRATCHDIR"

# Load required modules
module purge
module load uri/main
module load conda/latest

# Copy input files to scratch
cp /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.coding.gff3 $SCRATCHDIR
cp /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.transcripts.fasta $SCRATCHDIR
echo "[$(date)] Input files copied to scratch directory"
echo "[$(date)] Files in scratch directory:"
ls $SCRATCHDIR

# Set conda environment path
CONDA_ENV_PATH="$SCRATCHDIR/conda-envs/agat-env"

# Combine GFF3 with FASTA (include the ##FASTA tag)
cp Ahyacinthus.coding.gff3 tmp_ahya.gff3
echo "##FASTA" >> tmp_ahya.gff3
cat Ahyacinthus.transcripts.fasta >> tmp_ahya.gff3
mv tmp_ahya.gff3 combined_ahya.gff3
echo "[$(date)] Combined GFF3 with transcript FASTA into combined_ahya.gff3"

# Run bp_seqconvert with correct options
echo "[$(date)] Running bp_seqconvert conversion..."
conda run -p "$CONDA_ENV_PATH" \
bp_seqconvert --from gff3 --to genbank < combined_ahya.gff3 > ahya.gb
echo "[$(date)] GenBank file created."

# Validate output
if [ -f ahya.gb ]; then
    echo "[$(date)] Output file ahya.gb exists."
else
    echo "[$(date)] Error: Output file ahya.gb not created."
    exit 1
fi

# Copy result back to working directory
DEST_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann
cp ahya.gb "$DEST_DIR" && \
echo "[$(date)] Output file copied to $DEST_DIR" || \
{ echo "[$(date)] Error copying output file."; exit 1; }

# Final listing
echo "[$(date)] Files in destination:"
ls "$DEST_DIR"
