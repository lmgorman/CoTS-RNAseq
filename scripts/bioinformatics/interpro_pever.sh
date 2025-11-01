#Navigate to Porites spp. genome folder
cd /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por
#Download the P. evermanni protein fasta file
wget https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.pep.fa 
#Had to remove asterix from the downloaded file for interpro to work
see remove_asterix.sh script - https://github.com/lmgorman/CoTS-RNAseq/blob/main/scripts/bioinformatics/remove_asterix.sh
#Create output directory for interpro
mkdir /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro
mkdir /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output
cd /home/lucy_gorman_uri_edu/scripts
#Create following shell script in /home/lucy_gorman_uri_edu/scripts called interpro_pever.sh



#!/bin/bash
#SBATCH --job-name=annot-Pever-interpro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH -p uri-cpu
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -o annot-Pever-interpro-%j.out
#SBATCH -e annot-Pever-interpro-%j.error
# Load modules
echo "Loading programs" $(date)
module purge
module load uri/main
module load InterProScan/5.73-104.0-foss-2024a
# Use pre-assigned scratch directory
SCRATCHDIR=/scratch/workspace/lucy_gorman_uri_edu-lucyscratch
# Define input and output locations
INPUT_FILE=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/out.fa
BASENAME=$(basename "$INPUT_FILE")
FINAL_OUTPUT_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output
# Copy input file to scratch
cp "$INPUT_FILE" "$SCRATCHDIR/"
cd "$SCRATCHDIR"
# Run InterProScan
echo "Running InterProScan on $(hostname) at $(date)"
interproscan.sh \
    --input "$BASENAME" \
    --disable-precalc \
    --output-dir "$SCRATCHDIR" \
    --tempdir "$SCRATCHDIR/temp" \
    --cpu 10 \
    --goterms \
    --pathways
# Create final output directory if needed
mkdir -p "$FINAL_OUTPUT_DIR"
# Copy results back
echo "Copying results to $FINAL_OUTPUT_DIR"
cp "${BASENAME}"* "$FINAL_OUTPUT_DIR/" 2>/dev/null || echo "No direct-matching output, copying all results..."
cp "$SCRATCHDIR"/* "$FINAL_OUTPUT_DIR/"
# Optional: Keep scratch clean
echo "Cleaning up temporary files from scratch directory"
rm -rf "$SCRATCHDIR/temp"
echo "InterProScan completed at $(date)"
