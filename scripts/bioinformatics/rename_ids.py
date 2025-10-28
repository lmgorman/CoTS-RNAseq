#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import xml.etree.ElementTree as ET

def build_mapping(fasta_in, prefix="scaf"):
    """Create a mapping of old IDs to new short IDs based on FASTA"""
    mapping = {}
    for i, record in enumerate(SeqIO.parse(fasta_in, "fasta"), 1):
        new_id = f"{prefix}{i:07d}"
        mapping[record.id] = new_id
    return mapping

def rename_fasta(fasta_in, fasta_out, mapping):
    """Write renamed FASTA file"""
    with open(fasta_out, "w") as fout:
        for record in SeqIO.parse(fasta_in, "fasta"):
            record.id = mapping.get(record.id, record.id)
            record.description = ""
            SeqIO.write(record, fout, "fasta")

def rename_gff3(gff_in, gff_out, mapping):
    """Rename scaffolds in GFF3"""
    with open(gff_in) as fin, open(gff_out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            parts = line.strip().split("\t")
            parts[0] = mapping.get(parts[0], parts[0])
            fout.write("\t".join(parts) + "\n")

def rename_eggnog(egg_in, egg_out, mapping):
    """Rename gene/protein IDs in EggNOG mapper file"""
    with open(egg_in) as fin, open(egg_out, "w") as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            old_id = parts[0]
            # use mapping if available, otherwise keep original
            parts[0] = mapping.get(old_id, old_id)
            fout.write("\t".join(parts) + "\n")

def rename_interpro(interpro_in, interpro_out, mapping):
    """Rename protein IDs in InterProScan XML"""
    tree = ET.parse(interpro_in)
    root = tree.getroot()
    for protein in root.findall(".//protein"):
        pid = protein.attrib.get("id")
        if pid in mapping:
            protein.set("id", mapping[pid])
    tree.write(interpro_out)

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA, GFF3, EggNOG, and InterProScan IDs consistently")
    parser.add_argument("--fasta", required=True, help="Input FASTA")
    parser.add_argument("--gff", required=True, help="Input GFF3")
    parser.add_argument("--eggnog", required=True, help="Input EggNOG annotations")
    parser.add_argument("--interpro", required=True, help="Input InterProScan XML")
    parser.add_argument("--prefix", default="scaf", help="Prefix for new IDs")
    parser.add_argument("--outdir", default="renamed", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    fasta_out = outdir / (Path(args.fasta).stem + "_renamed.fa")
    gff_out = outdir / (Path(args.gff).stem + "_renamed.gff")
    egg_out = outdir / (Path(args.eggnog).stem + "_renamed.emapper.annotations")
    interpro_out = outdir / (Path(args.interpro).stem + "_renamed.xml")

    print("Building ID mapping from FASTA...")
    mapping = build_mapping(args.fasta, prefix=args.prefix)

    print(f"Renaming FASTA -> {fasta_out}")
    rename_fasta(args.fasta, fasta_out, mapping)

    print(f"Renaming GFF3 -> {gff_out}")
    rename_gff3(args.gff, gff_out, mapping)

    print(f"Renaming EggNOG -> {egg_out}")
    rename_eggnog(args.eggnog, egg_out, mapping)

    print(f"Renaming InterProScan XML -> {interpro_out}")
    rename_interpro(args.interpro, interpro_out, mapping)

    print(f"All files renamed consistently. Output directory: {outdir}")

if __name__ == "__main__":
    main()
