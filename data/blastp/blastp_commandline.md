Acropora hyacinthus blastp
`D:`
`cd "blast-2.17.0+\bin"`
makeblastdb -in D:\Acr_allexpressedgenes.fasta -dbtype prot -out D:\AllexpressedgenesAcro

Results
Building a new DB, current time: 12/23/2025 13:32:28
New DB name:   D:\AllexpressedgenesAcro
New DB title:  D:\Acr_allexpressedgenes.fasta
Sequence type: Protein
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 18258 sequences in 1.65395 seconds.


blastp -query "D:\toxprot_cnidarian_342.fasta" `
    -db "D:\AllexpressedgenesAcro" `
    -num_threads 4 `
    -evalue 1e-5 `
    -max_target_seqs 1 `
    -outfmt 5 `
    -out "D:\results.xml"

Results
87 hits, D:\results.xml

blastp -query "D:\toxprot_cnidarian_342.fasta" -db "D:\AllexpressedgenesAcro" -matrix BLOSUM62 -num_threads 4 -evalue 5e-2 -max_target_seqs 5000 -outfmt 5 -out "D:\Acr_results_relaxed.xml" -word_size 6 -gapopen 11 -gapextend 1 -seg no
Results
821 hits, D:\Acr_results_relaxed.xml
Strict results
blastp -query "D:\toxprot_cnidarian_342.fasta" -db "D:\AllexpressedgenesAcro" -matrix BLOSUM62 -num_threads 4 -evalue 1e-5 -max_target_seqs 5000 -outfmt 5 -out "D:\Acr_results_strict.xml" -word_size 6 -gapopen 11 -gapextend 1 -seg no -qcov_hsp_perc 70
348 results D:\Acr_results_strict.xml
Porites 
makeblastdb -in D:\Por_allexpressedgenes.fasta -dbtype prot -out D:\AllexpressedgenesPor
Relaxed parameters 
blastp -query "D:\toxprot_cnidarian_342.fasta" -db "D:\AllexpressedgenesPor" -matrix BLOSUM62 -num_threads 4 -evalue 5e-2 -max_target_seqs 5000 -outfmt 5 -out "D:\Por_results_relaxed.xml" -word_size 6 -gapopen 11 -gapextend 1 -seg no 
1350 results D:\Por_results_relaxed.xml
Stricter parameters 
blastp -query "D:\toxprot_cnidarian_342.fasta" -db "D:\AllexpressedgenesPor" -matrix BLOSUM62 -num_threads 4 -evalue 1e-5 -max_target_seqs 5000 -outfmt 5 -out "D:\Por_results_strict.xml" -word_size 6 -gapopen 11 -gapextend 1 -seg no -qcov_hsp_perc 70
