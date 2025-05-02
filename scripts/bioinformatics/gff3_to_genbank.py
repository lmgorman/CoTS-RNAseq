#Create conda environment in scratch directory in an interactive session
cd /scratch/workspace/lucy_gorman_uri_edu-lucyscratch
salloc -p cpu -c 2 --mem=5G --time 02:00:00
module load conda/latest
conda create --name Lucy python=3.7
conda activate Lucy





#Saved the following script in /scratch/workspace/lucy_gorman_uri_edu-lucyscratch/gff3_to_genbank.py
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import sys

# Input files
fasta_file = "Ahyacinthus.chrsV1.fasta"
gff_file = "Ahyacinthus.coding.gff3"
output_file = "ahya.gbk"

seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Required for GenBank output
for record in seq_dict.values():
    record.annotations["molecule_type"] = "DNA"

# Prepare output records
records = []

# Parse the GFF3 file manually
with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) != 9:
            continue  # skip malformed lines
        seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts

        # Parse attributes into a dictionary
        attr_dict = {}
        for attr in attributes.split(";"):
            if "=" in attr:
                key, value = attr.split("=", 1)
                attr_dict[key] = value

        # Build a SeqFeature
        location = FeatureLocation(int(start)-1, int(end), strand=1 if strand == "+" else -1)
        qualifiers = {k: [v] for k, v in attr_dict.items()}
        feat = SeqFeature(location=location, type=feature_type, qualifiers=qualifiers)

        # Add feature to the appropriate record
        if seq_id not in seq_dict:
            print(f"Warning: {seq_id} not found in FASTA. Skipping feature.")
            continue
        if not hasattr(seq_dict[seq_id], "features"):
            seq_dict[seq_id].features = []
        seq_dict[seq_id].features.append(feat)

# Now write GenBank-formatted records
records = list(seq_dict.values())
with open(output_file, "w") as out_handle:
    SeqIO.write(records, out_handle, "genbank")

print(f"✅ GenBank file written to: {output_file}")
