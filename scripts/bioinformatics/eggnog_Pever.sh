#!/bin/bash
#SBATCH --job-name=annot-Pever-eggnog
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH -p uri-cpu
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -o annot-Pever-eggnog-%j.out
#SBATCH -e annot-Pever-eggnog-%j.error
# Load modules
echo "Loading modules" $(date)
module purge
module load uri/main
module load all/eggnog-mapper/2.1.12-foss-2023a
# Use pre-assigned scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
# Define input and output
INPUT=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/out.fa #truncated names that match gff3 and fasta files
OUTPUT_BASENAME=pever_eggnog_trunc
FINAL_OUTPUT_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/eggnog
# Copy input file to scratch
cp "$INPUT" "$SCRATCHDIR/"
cd "$SCRATCHDIR"
# Run eggnog-mapper
echo "Running eggnog-mapper on $(hostname) at $(date)"
/modules/uri_apps/software/eggnog-mapper/2.1.12-foss-2023a/bin/emapper.py \
  -i $(basename "$INPUT") \
  --itype proteins \
  -m diamond \
  --output "$OUTPUT_BASENAME" \
  --cpu 10 \
  --go_evidence non-electronic \
  --data_dir /modules/uri_apps/software/eggnog-mapper/2.1.12-foss-2023a/lib/python3.11/site-packages/data \
  --override
# Copy results back
echo "Copying results back to $FINAL_OUTPUT_DIR"
mkdir -p "$FINAL_OUTPUT_DIR"
cp "${OUTPUT_BASENAME}"* "$FINAL_OUTPUT_DIR/"
# Optional: Keep scratch clean
echo "Cleaning up scratch"
rm -rf "$SCRATCHDIR"
echo "Job complete at $(date)"
