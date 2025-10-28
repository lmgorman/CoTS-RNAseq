#!/bin/bash
#SBATCH --job-name=agat-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o agat-pever-%j.out
#SBATCH -e agat-pever-%j.error

set -e  # Exit on first error

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR

echo "[$(date)] Job started in $SCRATCHDIR"

# Load required modules
module purge
module load uri/main
module load conda/latest

# Copy input files to scratch
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/Porites_evermanni_v1_FIXED.gff $SCRATCHDIR
echo "[$(date)] Input files copied to scratch directory"
echo "[$(date)] Files in scratch directory:"
ls $SCRATCHDIR

# Set conda environment path
CONDA_ENV_PATH="$SCRATCHDIR/conda-envs/agat-env"

# Run bp_seqconvert with correct options
echo "[$(date)] Running bp_seqconvert conversion..."
conda run -p "$CONDA_ENV_PATH" \
bp_seqconvert --from gff3 --to genbank < Porites_evermanni_v1_FIXED.gff > Pever.gb
echo "[$(date)] GenBank file created."

# Validate output
if [ -f Pever.gb ]; then
    echo "[$(date)] Output file Pever.gb exists."
else
    echo "[$(date)] Error: Output file Pever.gb not created."
    exit 1
fi

# Copy result back to working directory
DEST_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/
cp ahya.gb "$DEST_DIR" && \
echo "[$(date)] Output file copied to $DEST_DIR" || \
{ echo "[$(date)] Error copying output file."; exit 1; }

# Final listing
echo "[$(date)] Files in destination:"
ls "$DEST_DIR"
