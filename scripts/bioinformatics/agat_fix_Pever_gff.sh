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

# Load required modules
module purge
module load uri/main
module load conda/latest

# -------------------------
# Activate conda env from working directory
# -------------------------
WORKDIR="/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever"
ENV_PATH="$WORKDIR/envs/agat-env"

# activate your local env
source activate "$ENV_PATH"

cd "$WORKDIR"

echo "Running AGAT from $(pwd) at $(date)"
echo "Using conda environment at: $ENV_PATH"


# Run AGAT
# -------------------------
echo "Starting AGAT GFF repair at $(date)"

# If you have a big GFF file, add --log to capture AGAT messages
agat_convert_sp_gff2gff.pl \
    --gff Porites_evermanni_v1_CORRECT.gff \
    -o Porites_evermanni_v1_FIXED.gff \
    --log agat_fix.log

echo "AGAT finished at $(date)"
