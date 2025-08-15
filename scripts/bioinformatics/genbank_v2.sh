#!/bin/bash
#SBATCH --job-name=Ahya_gbff_to_gff_gtf
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH -p uri-cpu
#SBATCH --time=24:00:00
#SBATCH -o Ahya_gbff_to_gff_gtf-%j.out
#SBATCH -e Ahya_gbff_to_gff_gtf-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1

# Load modules
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0
module load gffread/0.12.7

#Convert to gff3
gffread /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff -o GCA_964291705.1_jaAcrHyac4.1_genomic_2.gff
#convert to gtf
gffread /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff -T -o GCA_964291705.1_jaAcrHyac4.1_genomic_2.gtf
