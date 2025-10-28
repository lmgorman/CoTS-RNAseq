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

set -e  # Exit immediately if any command fails

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# Load modules
module purge
module load uri/main
module load funannotate/1.8.17

# Define Apptainer container
FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"

# Copy input data to scratch
echo "[$(date)] Copying input data to scratch..."
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/Porites_evermanni_v1.fa $SCRATCHDIR
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/Porites_evermanni_v1_FIXED.gff $SCRATCHDIR
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output/Porites_evermanni_v1_clean.annot.pep.fa.xml $SCRATCHDIR/iprscan.xml
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/eggnog/pever_eggnog.emapper.annotations $SCRATCHDIR/eggnog.annotations

# Run Funannotate annotate
echo "[$(date)] Starting Funannotate annotation..."
apptainer run "$FUNANNOTATE_SIF" funannotate annotate \
  --gff $SCRATCHDIR/Porites_evermanni_v1_FIXED.gff \
  --fasta $SCRATCHDIR/Porites_evermanni_v1.fa \
  -s "Porites evermanni" \
  -o $SCRATCHDIR/output \
  --iprscan $SCRATCHDIR/iprscan.xml \
  --eggnog $SCRATCHDIR/eggnog.annotations \
  --busco_db metazoa \
  --force \
  --cpus 10

echo "[$(date)] Renaming long scaffold IDs in annotation outputs..."
for f in $SCRATCHDIR/output/*.gff3 $SCRATCHDIR/output/*.faa $SCRATCHDIR/output/*.fna $SCRATCHDIR/output/*.annotation.txt; do
    python /path/to/rename_ids_generic.py "$f" "${f%.txt}_shortIDs.txt"
done
echo "[$(date)] Scaffold renaming complete."

# Patch LOCUS lines in GenBank files immediately after Funannotate finishes
echo "[$(date)] Patching LOCUS lines in GenBank files..."
for gbk in $SCRATCHDIR/output/*.gbk; do
    sed -i 's/ bp   DNA/ 1000 bp   DNA/' "$gbk"
done

# Copy results back to work
echo "[$(date)] Copying results back to /work..."
rsync -av "$SCRATCHDIR/output/" "/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/funannotate/"

echo "[$(date)] Annotation complete. Job finished successfully."
