#!/bin/bash
#SBATCH --job-name=funannotate-ahya
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH -p uri-cpu
#SBATCH -o funannotate-ahya-%j.out
#SBATCH -e funannotate-ahya-%j.error

# Load Funannotate module or activate conda environment
module purge
module load all/funannotate/1.8.15-foss-2023a
# or, if using conda: source activate funannotate-env

# Run annotation
funannotate annotate \
  --genbank /work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.gb \
  -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahyacinthus_funannotate \
  --iprscan /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/interpro/output/Ahyacinthus.proteins.fasta_1.xml \
  --eggnog /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/eggnog/X.emapper.annotations \
  --busco_db metazoa \
  --cpus 10
