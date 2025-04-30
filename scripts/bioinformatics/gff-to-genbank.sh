#!/bin/bash
#SBATCH --job-name=gff-to-genbank
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o gff2gb-%j.out
#SBATCH -e gff2gb-%j.error

set -e  # Exit on error

# Define directories
SCRATCHDIR=/scratch/workspace/lucy_gorman_uri_edu-lucyscratch
WORKDIR=/work/pi_hputnam_uri_edu/refs/Ahyacinthus_genome/Ahyacinthus_genome_V1
DEST_DIR=/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann

# Go to scratch
cd $SCRATCHDIR
echo "[$(date)] Job started in $SCRATCHDIR"

# Load Python module
module purge
module load uri/main
module load python/3.12.3
python3.12 -m pip install --user biopython bcbio-gff
python3.12 -c "from BCBio import GFF; from Bio import SeqIO; print('All good!')"

# Copy input files
cp $WORKDIR/Ahyacinthus.coding.gff3 .
cp $WORKDIR/Ahyacinthus.chrsV1.fasta .
echo "[$(date)] Files copied:"
ls -lh

# Run Biopython script
echo "[$(date)] Running Biopython GFF to GenBank conversion..."
python3.12 <<EOF
from BCBio import GFF
from Bio import SeqIO

in_gff = "Ahyacinthus.coding.gff3"
in_fasta = "Ahyacinthus.chrsV1.fasta"
out_gbk = "ahya.gbk"

fasta_dict = SeqIO.to_dict(SeqIO.parse(open(in_fasta), "fasta"))

with open(in_gff) as gff_handle, open(out_gbk, "w") as gbk_handle:
    GFF.write(GFF.parse(gff_handle, base_dict=fasta_dict), gbk_handle)
EOF

# Confirm output
if [ -f ahya.gbk ]; then
    echo "[$(date)] ahya.gbk successfully created."
else
    echo "[$(date)] ERROR: ahya.gbk was not created."
    exit 1
fi

# Copy result to destination
cp ahya.gbk $DEST_DIR && \
echo "[$(date)] File copied to $DEST_DIR" || \
{ echo "[$(date)] Copy failed."; exit 1; }

# Final listing
echo "[$(date)] Files in destination:"
ls -lh $DEST_DIR

