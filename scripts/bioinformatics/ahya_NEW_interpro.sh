#!/bin/bash
#SBATCH --job-name=annot-AhyaNEW-interpro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH -p uri-cpu
#SBATCH --time=96:00:00
#SBATCH --constraint=avx512
#SBATCH -o annot-AhyaNEW-interpro-%j.out
#SBATCH -e annot-AhyaNEW-interpro-%j.error

# Load modules
echo "Loading programs" $(date)
module purge
module load uri/main
module load InterProScan/5.73-104.0-foss-2024a

# Scratch and paths
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
# Need protein file to input
INPUT_FILE=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.fna
BASENAME=$(basename "$INPUT_FILE")
FINAL_OUTPUT_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/interpro

# Copy input to scratch
cp "$INPUT_FILE" "$SCRATCHDIR/"
cd "$SCRATCHDIR"

# Run InterProScan
echo "Running InterProScan on $(hostname) at $(date)"
interproscan.sh \
    --input "$BASENAME" \
    --type n \
    --disable-precalc \
    --output-dir "$SCRATCHDIR" \
    --tempdir "$SCRATCHDIR/temp" \
    --cpu 10 \
    --formats tsv \
    --goterms \
    --pathways

# Create output dir if needed
mkdir -p "$FINAL_OUTPUT_DIR"

# Copy result files (safe and filtered)
echo "Copying results to $FINAL_OUTPUT_DIR"
cp "$SCRATCHDIR"/*.tsv "$FINAL_OUTPUT_DIR/" 2>/dev/null || echo "No TSV output found."

# Clean up
echo "Cleaning up temporary files from scratch directory"
rm -rf "$SCRATCHDIR/temp"

echo "InterProScan completed at $(date)"
