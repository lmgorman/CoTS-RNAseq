#!/bin/bash
#SBATCH --job-name=agat_fix_Pever_gff
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o agat_fix_Pever_gff-%j.out
#SBATCH -e agat_fix_Pever_gff-%j.error

# -------------------------
# Define scratch directory
# -------------------------
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# -------------------------
# Load required modules
# -------------------------
module purge
module load uri/main
module load conda/latest


# Ensure conda functions are available
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/conda-envs/agat-env


# -------------------------
# Copy input files to scratch
# -------------------------
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/Porites_evermanni_v1_CORRECT.gff $SCRATCHDIR
echo "[$(date)] Input files copied to scratch directory"
echo "[$(date)] Files in scratch directory:"
ls $SCRATCHDIR

# -------------------------
# Run AGAT
# -------------------------
OUTPUT_DIR="$SCRATCHDIR/AGAT_fixed_output"
mkdir -p "$OUTPUT_DIR"

agat_convert_sp_gff2gff.pl \
  --gff Porites_evermanni_v1_CORRECT.gff \
  -o "$OUTPUT_DIR/Porites_evermanni_v1_FIXED.gff" \
  --log "$OUTPUT_DIR/agat_fix.log"

echo "[$(date)] Finished AGAT run. Output in $OUTPUT_DIR"

# -------------------------
# Optional: move results back to permanent storage
# -------------------------
# cp -r "$OUTPUT_DIR" /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/
