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

# Define the Apptainer container path
FUNANNOTATE_SIF="/modules/opt/linux-ubuntu24.04-x86_64/funannotate/1.8.17/funannotate-1.8.17.sif"

echo "[$(date)] Loading funannotate module..."

# Load Funannotate module or activate conda environment
module purge
module load uri/main
module load funannotate/1.8.17
# or, if using conda: source activate funannotate-env

# Check if the module was successfully loaded
echo "[$(date)] Loaded modules: $(module list)"

# Run annotation
apptainer run "$FUNANNOTATE_SIF" funannotate annotate \
  -i /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_funann/ \
  -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahyaannotate_results \
  --iprscan /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/interpro/output/Ahyacinthus.proteins.fasta_1.xml \
  --eggnog /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/eggnog/ahya_eggnog.emapper.annotations \
  --busco_db metazoa \
  --cpus 10
echo "[$(date)] Annotation complete. Job finished successfully."
