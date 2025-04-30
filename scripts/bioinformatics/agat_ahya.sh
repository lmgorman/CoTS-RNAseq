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

set -e  # Stop on first error

# Define scratch directory
SCRATCHDIR=/scratch/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR

echo "[$(date)] Job started in $SCRATCHDIR"

# Load required modules
module purge
module load uri/main
module load conda/latest

# Copy files to scratch directory, checking if they exist first
SOURCE_GFF="/work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.coding.gff3"
SOURCE_FASTA="/work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.proteins.fasta"

if [[ -f "$SOURCE_GFF" && -f "$SOURCE_FASTA" ]]; then
    cp $SOURCE_GFF $SCRATCHDIR
    cp $SOURCE_FASTA $SCRATCHDIR
    echo "[$(date)] Input files copied to scratch directory"
else
    echo "[$(date)] Error: One or both input files do not exist. Exiting."
    exit 1
fi

echo "[$(date)] Files in scratch directory: $(ls $SCRATCHDIR)"

# Initialize conda
source $(conda info --base)/etc/profile.d/conda.sh

# Create environment if it doesn't exist (check for AGAT availability first)
CONDA_ENV_PATH="$SCRATCHDIR/conda-envs/agat-env"
if ! conda list | grep -q agat; then
    echo "[$(date)] Creating AGAT conda environment..."
    conda create -y -p "$CONDA_ENV_PATH" -c bioconda -c conda-forge agat
else
    echo "[$(date)] AGAT environment already exists."
fi

# Activate it properly
conda activate $SCRATCHDIR/conda-envs/agat-env
echo "[$(date)] AGAT environment activated."

# List installed packages to check if AGAT is installed
echo "[$(date)] Listing installed packages in the environment:"
conda list

# Check if AGAT tool exists in the environment
echo "[$(date)] AGAT location: $(which agat_convert_sp_gff2gbk.pl)"

if which agat_convert_sp_gff2gbk.pl; then
    echo "[$(date)] Running AGAT conversion..."
    $SCRATCHDIR/conda-envs/agat-env/bin/agat_convert_sp_gff2gbk.pl \
      --gff Ahyacinthus.coding.gff3 \
      --fasta Ahyacinthus.proteins.fasta \
      --output ahya.gb
    echo "[$(date)] GenBank file successfully created"
else
    echo "[$(date)] Error: AGAT tool not found. Exiting."
    exit 1
fi

# Check if the GenBank file was created
if [ -f ahya.gb ]; then
    echo "[$(date)] GenBank file created successfully."
else
    echo "[$(date)] Error: GenBank file not created."
    exit 1
fi

# Copy output GenBank file back to working directory
DEST_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann
if cp ahya.gb $DEST_DIR; then
    echo "[$(date)] Output file copied to $DEST_DIR"
else
    echo "[$(date)] Error: Could not copy GenBank file to $DEST_DIR"
    exit 1
fi

# Final directory listing for confirmation
echo "[$(date)] Files in destination directory: $(ls $DEST_DIR)"
