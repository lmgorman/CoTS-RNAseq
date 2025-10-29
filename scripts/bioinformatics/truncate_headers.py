#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 5:
    print("Usage: python truncate_headers.py input.fasta input.gff3 output.fasta output.gff3")
    sys.exit(1)

fasta_in, gff_in, fasta_out, gff_out = sys.argv[1:]

# -------- STEP 1: Parse FASTA and rename headers --------
id_map = {}

with open(fasta_in) as fin, open(fasta_out, "w") as fout:
    counter = 1
    for line in fin:
        if line.startswith(">"):
            orig_id = line[1:].strip().split()[0]
            new_id = f"scaffold_{counter}"
            id_map[orig_id] = new_id
            fout.write(f">{new_id}\n")
            counter += 1
        else:
            fout.write(line)

# -------- STEP 2: Rewrite GFF3 using new names --------
def replace_seqid(line):
    if line.startswith("#") or line.strip() == "":
        return line
    parts = line.split("\t")
    if len(parts) < 9:
        return line
    old_seqid = parts[0]
    if old_seqid in id_map:
        parts[0] = id_map[old_seqid]
    return "\t".join(parts)

with open(gff_in) as fin, open(gff_out, "w") as fout:
    for line in fin:
        fout.write(replace_seqid(line))

# -------- STEP 3: Write mapping file (for reference) --------
with open("seqid_mapping.tsv", "w") as mapout:
    mapout.write("Original_ID\tNew_ID\n")
    for k, v in id_map.items():
        mapout.write(f"{k}\t{v}\n")

print(f"✅ Done!\nRenamed {len(id_map)} sequences.")
print("Mapping written to seqid_mapping.tsv")
