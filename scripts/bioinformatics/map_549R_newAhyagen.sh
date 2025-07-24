#!/bin/bash
#SBATCH --job-name=Ahya_align_star
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH -D work/pi_hputnam_uri_edu/20250107_COTS_LG/

#load modules
echo "Loading programs" $(date)
module load uri/main
module load STAR/2.7.11b-GCC-12.3.0

echo "Starting read alignment." $(date)

STAR --runMode alignReads \
--genomeDir  /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1 \
--runThreadN 10 \
--readFilesCommand zcat \
--readFilesIn /work/pi_hputnam_uri_edu/20250107_COTS_LG/fastp_trimmed/Ahya/549R_R1_001_trimmed.fastq.gz /work/pi_hputnam_uri_edu/20250107_COTS_LG/fastp_trimmed/Ahya/549R_R2_001_trimmed.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFileNamePrefix ${i%_newgenome.fastq.gz}
done

echo "Alignment of Trimmed Seq data complete." $(date)
