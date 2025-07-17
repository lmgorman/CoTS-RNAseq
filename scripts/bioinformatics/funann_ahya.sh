#!/bin/bash
#SBATCH --job-name=funannotate-ahya
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o funannotate-ahya-%j.out
#SBATCH -e funannotate-ahya-%j.error

set -e  # Exit immediately if a command exits with a non-zero status

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# Load necessary modules
module purge
module load uri/main
module load funannotate/1.8.17

# Define Apptainer container
FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"

# Copy input data to scratch
echo "[$(date)] Copying input data to scratch..."
cp -r /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_funann $SCRATCHDIR/input

# Copy interpro and eggnog files to scratch (optional but recommended for full scratch use)
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/interpro/output/Ahyacinthus.proteins.fasta_1.xml $SCRATCHDIR/iprscan.xml
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/eggnog/ahya_eggnog.emapper.annotations $SCRATCHDIR/eggnog.annotations

# Run funannotate inside Apptainer from scratch
echo "[$(date)] Starting funannotate..."
apptainer run "$FUNANNOTATE_SIF" funannotate annotate \
  -i $SCRATCHDIR/input \
  -o $SCRATCHDIR/output \
  --iprscan $SCRATCHDIR/iprscan.xml \
  --eggnog $SCRATCHDIR/eggnog.annotations \
  --busco_db metazoa \
  --cpus 10

# Copy results back to work
echo "[$(date)] Copying results back to /work..."
cp -r $SCRATCHDIR/output /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahyaannotate_results

echo "[$(date)] Annotation complete. Job finished successfully."
