#!/bin/bash
#SBATCH --job-name=funannotate-predict-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o funannotate-predict-pever-%j.out
#SBATCH -e funannotate-predict-pever-%j.error

set -e  # Exit immediately if any command fails

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# Load modules
module purge
module load uri/main
module load funannotate/1.8.17

# Check if the module was successfully loaded
echo "[$(date)] Loaded modules: $(module list)"

# Turn genome into softmasked repeat genome
apptainer run "$FUNANNOTATE_SIF" funannotate mask \
             -i /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta \
             -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_funann/Ahyacinthus_sm.chrsV1.fasta 
             
# Download BUSCO database if it does not exist

# Check if the BUSCO database is already downloaded
if [ ! -d "$FUNANNOTATE_DB/metazoa" ]; then
    apptainer run "$FUNANNOTATE_SIF" funannotate setup -b metazoa -d "$FUNANNOTATE_DB"
else
    echo "[$(date)] metazoa BUSCO database already exists."
fi

# Define Apptainer container and database
FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"
FUNANNOTATE_DB=/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/database

# Copy genome and optional RNA evidence to scratch
echo "[$(date)] Copying input data to scratch..."
cp /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/Porites_evermanni_v1.fa $SCRATCHDIR

# Create output directory
OUTDIR=$SCRATCHDIR/predict_output
rm -rf $OUTDIR
mkdir -p $OUTDIR

# Run Funannotate predict
echo "[$(date)] Starting Funannotate predict..."
apptainer run "$FUNANNOTATE_SIF" funannotate predict \
  -i $SCRATCHDIR/Porites_evermanni_v1.fa \
  -o $OUTDIR \
  -s "Porites evermanni" \
  --busco_db metazoa \
  --cpus 16 \
  --optimize_augustus \
  --strain v1 \
  --force \
  --tmpdir $SCRATCHDIR/tmp

echo "[$(date)] Predict finished. Copying results back..."
rsync -av "$OUTDIR/" "/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/funannotate_predict/"

echo "[$(date)] Funannotate predict completed successfully."
