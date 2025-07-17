#!/bin/bash
#SBATCH --job-name=funannotate-ahya
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o funpred-opti-ahya-%j.out
#SBATCH -e funpred-opti-ahya-%j.error

set -e  # Exit immediately if a command exits with a non-zero status

# Define scratch directory
SCRATCHDIR=/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# Define the Apptainer container path
FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"

# Define the FUNANNOTATE_DB path to the newly downloaded database
export FUNANNOTATE_DB="/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/funannotate_databases"
echo "[$(date)] Loading funannotate module..."

# Load Funannotate module or activate conda environment
module purge
module load uri/main
module load funannotate/1.8.17

# Check if the module was successfully loaded
echo "[$(date)] Loaded modules: $(module list)"
          
# Download BUSCO database if it does not exist

# Check if the BUSCO database is already downloaded
if [ ! -d "$FUNANNOTATE_DB/metazoa" ]; then
    apptainer run "$FUNANNOTATE_SIF" funannotate setup -b metazoa -d "$FUNANNOTATE_DB"
else
    echo "[$(date)] metazoa BUSCO database already exists."
fi

# Run annotation
apptainer run "$FUNANNOTATE_SIF" funannotate predict \
            -i /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_funann/Ahyacinthus_sm.chrsV1.fasta  \
            --species "Acropora hyacinthus" \
            -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_pred_opti \
            --protein_evidence /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.proteins.fasta \
            --transcript_evidence /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.coding.gff3 \
            --busco_db "$FUNANNOTATE_DB/metazoa" \
            --optimize_augustus
            --force            
            --cpus 8
            --verbose
            
echo "[$(date)] Prediction complete. Job finished successfully."
