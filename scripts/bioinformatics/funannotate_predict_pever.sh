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

#----------------------------------------------------------
# 1. Setup environment
#----------------------------------------------------------
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

module purge
module load uri/main
module load funannotate/1.8.17

FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"
FUNANNOTATE_DB="$SCRATCHDIR/funannotate_databases"
export FUNANNOTATE_DB

echo "[$(date)] Using database path: $FUNANNOTATE_DB"

#----------------------------------------------------------
# 2. Mask repeats in genome
#----------------------------------------------------------
GENOME_ORIG="/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por-ever/Porites_evermanni_v1.fa"
GENOME_MASKED="$SCRATCHDIR/Porites_evermanni_v1_sm.fa"

echo "[$(date)] Running funannotate mask..."
apptainer run "$FUNANNOTATE_SIF" funannotate mask \
  -i "$GENOME_ORIG" \
  -o "$GENOME_MASKED"

#----------------------------------------------------------
# 3. Ensure BUSCO DB exists
#----------------------------------------------------------
if [ ! -d "$FUNANNOTATE_DB/metazoa" ]; then
  echo "[$(date)] Downloading BUSCO metazoa database..."
  apptainer run "$FUNANNOTATE_SIF" funannotate setup -b metazoa -d "$FUNANNOTATE_DB"
else
  echo "[$(date)] metazoa BUSCO database already exists."
fi

#----------------------------------------------------------
# 4. Run Funannotate predict
#----------------------------------------------------------
OUTDIR=$SCRATCHDIR/predict_output
rm -rf $OUTDIR
mkdir -p $OUTDIR

echo "[$(date)] Starting Funannotate predict..."
apptainer run "$FUNANNOTATE_SIF" funannotate predict \
  -i "$GENOME_MASKED" \
  -o "$OUTDIR" \
  -s "Porites evermanni" \
  --busco_db metazoa \
  --cpus 16 \
  --optimize_augustus \
  --strain v1 \
  --force \
  --tmpdir "$SCRATCHDIR/tmp"

#----------------------------------------------------------
# 5. Copy results back
#----------------------------------------------------------
echo "[$(date)] Predict finished. Copying results back..."
rsync -av "$OUTDIR/" "/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/funannotate_predict/"

echo "[$(date)] Funannotate predict completed successfully."
