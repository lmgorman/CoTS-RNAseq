#!/bin/bash
#SBATCH --job-name=clean_fasta
#SBATCH --output=clean_fasta.out
#SBATCH --error=clean_fasta.err
#SBATCH --time=01:00:00           
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# Load modules if needed
# module load uri/main

# Run the command to remove '*' characters
sed 's/\*//g' /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/Porites_evermanni_v1.annot.pep.fa \
  > /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/Porites_evermanni_v1_clean.annot.pep.fa

#############
Submitted above script on 25.04.25 under job id 33156060
