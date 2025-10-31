import pandas as pd

# Load EggNOG-mapper TSV
emapper_file = "eggnog.emapper.annotations"
df = pd.read_csv(emapper_file, sep="\t", comment="#")

rows = []

for idx, row in df.iterrows():
    gene = row['#query']
    
    if pd.notna(row['GOs']):
        for go in row['GOs'].split(','):
            rows.append([gene, 'GO', go])
    
    if pd.notna(row['EC']):
        for ec in row['EC'].split(','):
            rows.append([gene, 'EC', ec])
    
    if pd.notna(row['KEGG_ko']):
        for ko in row['KEGG_ko'].split(','):
            rows.append([gene, 'KEGG', ko])
    
    if pd.notna(row['PFAMs']):
        for pfam in row['PFAMs'].split(','):
            rows.append([gene, 'PFAM', pfam])

# Save Funannotate-compatible TSV
funannotate_tsv = "funannotate_annotations.tsv"
pd.DataFrame(rows, columns=['gene_id', 'source', 'term']).to_csv(funannotate_tsv, sep="\t", index=False, header=False)
