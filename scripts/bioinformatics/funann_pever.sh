#!/bin/bash
#SBATCH --job-name=funannotate-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o funannotate-pever-%j.out
#SBATCH -e funannotate-pever-%j.error

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
cp -r /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/Porites_evermanni_v1.fa $SCRATCHDIR
cp -r /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/Porites_evermanni_v1_FIXED.gff $SCRATCHDIR
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output/Porites_evermanni_v1_clean.annot.pep.fa.xml $SCRATCHDIR/iprscan.xml
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/eggnog/pever_eggnog.emapper.annotations $SCRATCHDIR/eggnog.annotations

echo "[$(date)] Shortening FASTA headers in original genome..."
awk '/^>/ {printf(">scaf%07d\n", ++i); next} {print}' $SCRATCHDIR/Porites_evermanni_v1.fa > $SCRATCHDIR/Porites_evermanni_v1.short.fa
mv $SCRATCHDIR/Porites_evermanni_v1.short.fa $SCRATCHDIR/Porites_evermanni_v1.fa

# Run funannotate inside Apptainer from scratch
echo "[$(date)] Starting funannotate..."
apptainer run "$FUNANNOTATE_SIF" funannotate annotate \
  --gff $SCRATCHDIR/Porites_evermanni_v1_FIXED.gff \
  --fasta $SCRATCHDIR/Porites_evermanni_v1.fa \
  -s "Porites evermanni" \
  -o $SCRATCHDIR/output \
  --iprscan $SCRATCHDIR/iprscan.xml \
  --eggnog $SCRATCHDIR/eggnog.annotations \
  --busco_db metazoa \
  --cpus 10

# FIX: Patch LOCUS lines in .gbk so Biopython can parse it
echo "[$(date)] Patching LOCUS lines in GenBank file..."
sed -i 's/ bp   DNA/ 1000 bp   DNA/' $SCRATCHDIR/output/*.gbk

echo "[$(date)] Shortening scaffold headers to <=16 characters..."
sed -i -E 's/(LOCUS       .{16}).*/\1/' $SCRATCHDIR/output/*.gbk

# Re-run only to regenerate Basename.annotations.txt (cached data, runs fast)
echo "[$(date)] Re-running funannotate annotate to regenerate annotations.txt..."
apptainer run "$FUNANNOTATE_SIF" funannotate annotate \
  --input $SCRATCHDIR/output \
  --force \
  --cpus 10

# Copy results back to work
echo "[$(date)] Copying results back to /work..."
rsync -av "$SCRATCHDIR/output/" "/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/funannotate/"

echo "[$(date)] Annotation complete. Job finished successfully."
