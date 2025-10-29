#!/bin/bash
#SBATCH --job-name=funannotate-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o funannotate-pever-trunc-%j.out
#SBATCH -e funannotate-pever-trunc-%j.error

set -e  # Exit immediately if any command fails

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch

# Clean output completely and create it
mkdir -p $SCRATCHDIR
rm -rf $SCRATCHDIR/output
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
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output/Porites_evermanni_v1_clean.annot.pep.fa.xml $SCRATCHDIR/iprscan.xml
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/eggnog/pever_eggnog.emapper.annotations $SCRATCHDIR/eggnog.annotations

# Ensure FASTA and GFF3 exist in scratch
if [[ ! -f $SCRATCHDIR/truncated_Porites_evermanni_v1.fasta ]] || [[ ! -f $SCRATCHDIR/truncated_Porites_evermanni_v1_FIXED.gff3 ]]; then
    echo "Error: FASTA or GFF3 not found in scratch."
    exit 1
fi

echo "[$(date)] Starting Funannotate annotation..."
apptainer run --bind $SCRATCHDIR:$SCRATCHDIR "$FUNANNOTATE_SIF" funannotate annotate \
  --gff $SCRATCHDIR/truncated_Porites_evermanni_v1_FIXED.gff3 \
  --fasta $SCRATCHDIR/truncated_Porites_evermanni_v1.fasta \
  -o $SCRATCHDIR/output \
  --iprscan $SCRATCHDIR/iprscan.xml \
  --eggnog $SCRATCHDIR/eggnog.annotations \
  --busco_db metazoa \
  --header_length 50 \
  --force \
  --cpus 10 

# Copy results back to work
echo "[$(date)] Copying results back to /work..."
rsync -av "$SCRATCHDIR/output/" "/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/funannotate/"

echo "[$(date)] Annotation complete. Job finished successfully."
