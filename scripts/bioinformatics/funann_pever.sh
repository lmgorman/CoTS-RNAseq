#!/bin/bash
#SBATCH --job-name=funannotate-pever
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH -p uri-cpu
#SBATCH -o funannotate-pever-%j.out
#SBATCH -e funannotate-pever-%j.error

# Load Funannotate module or activate conda environment
module purge
module load uri/main
module load all/funannotate/1.8.17
# or, if using conda: source activate funannotate-env

#funannotate import \
  --gff /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/Porites_evermanni_v1.annot.gff \
  --fasta /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/Porites_evermanni_v1.annot.pep.fa \
  --genbank pever.gb \
  --species "Porites evermanni" \
  --cpus 10


# Run annotation
funannotate annotate \
  --genbank /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/pever.gb \
  -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/por_ann/pever_funannotate \
  --iprscan /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/interpro/output/Porites_evermanni_v1_clean.annot.pep.fa.xml\
  --eggnog /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por/eggnog/pever_eggnog.emapper.annotations \
  --busco_db metazoa \
  --cpus 10
