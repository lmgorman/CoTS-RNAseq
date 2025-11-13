# 12.11.25
Created new directory under
`/work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn`


Created new genome directory to save the genomes into
`/work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome`

To use the method from
(https://github.com/pedronachtigall/ToxCodAn-Genome)

salloc -p uri-cpu -c 4 --mem=32G
module load uri/main
module load conda/latest

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home/lucy_gorman_uri_edu/.conda/envs/ToxcodanGenome



Downloaded their anthozoan CDS
https://github.com/pedronachtigall/ToxCodAn-Genome/blob/main/Databases/Anthozoa_db.fasta


# 13.11.25
## P. aus
Ran P.aus genome on default settings

    python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10
    
  ToxCodAn-Genome stopped to run, due to the following issue:
  No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

  ## P.rus
  python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Prus_genome.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

  ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

## P. compressa
  Porites_compressa.fasta
   python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Porites_compressa.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

   ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.
    
## P. lutea 
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/plut_final_2.1.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

 ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.
    
## P. astreoides
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Pastreoides.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

 ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.
    

## A. hyacinthus 
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Ahyacinthus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

2025-11-13 14:30:32 >>>> searching for toxin regions in genome using the database file...
2025-11-13 14:30:38 >>>> retrieving and filtering matched regions...
2025-11-13 14:30:38 >>>> processing overlapped regions...
2025-11-13 14:30:38 >>>> generating final output...
FASTA index file /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Ahyacinthus_genomic.fna.fai created.
        >>> Number of toxin loci identified in the genome:
                ACRMI -> 1
                TOTAL -> 1
2025-11-13 14:30:40 >>>> comparing annotated toxins to the uniprot/toxprot database...
        Uniprot/Toxprot database -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta
        >>> Check the uniprot/toxprot report: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/uniprot_report.txt
2025-11-13 14:30:40 >>>> ToxCodAn-Genome finished!
        >>> Check the final annotation output: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/toxin_annotation.gtf

## A. millepora 
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Amil_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Amil -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10


        Genome file -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Amil_genomic.fna
        Database file -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta
        Output folder -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Amil/
        Minimum percent identity -> 20
        Minimum gene size -> 20
        Maximum gene size -> 50000
        Minimum size of a CDS -> 20
        Number of threads -> 10
2025-11-13 14:32:56 >>>> searching for toxin regions in genome using the database file...
2025-11-13 14:33:02 >>>> retrieving and filtering matched regions...
2025-11-13 14:33:02 >>>> processing overlapped regions...
2025-11-13 14:33:02 >>>> generating final output...
FASTA index file /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Amil_genomic.fna.fai created.
        >>> Number of toxin loci identified in the genome:
                ACRMI -> 1
                TOTAL -> 1
2025-11-13 14:33:04 >>>> comparing annotated toxins to the uniprot/toxprot database...
        Uniprot/Toxprot database -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta
        >>> Check the uniprot/toxprot report: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Amil/uniprot_report.txt
2025-11-13 14:33:04 >>>> ToxCodAn-Genome finished!
        >>> Check the final annotation output: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Amil/toxin_annotation.gtf

## A. cervicornis 
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Acervicornis_genome.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Acerv -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10


ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

 ## Adigitifera
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Adigitifera_genome.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Adig -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

2025-11-13 14:36:42 >>>> starting ToxCodAn-Genome (v1.0 January 2023)...
        Genome file -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Adigitifera_genome.fna
        Database file -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta
        Output folder -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Adig/
        Minimum percent identity -> 20
        Minimum gene size -> 20
        Maximum gene size -> 50000
        Minimum size of a CDS -> 20
        Number of threads -> 10
2025-11-13 14:36:42 >>>> searching for toxin regions in genome using the database file...
2025-11-13 14:36:47 >>>> retrieving and filtering matched regions...
2025-11-13 14:36:47 >>>> processing overlapped regions...
2025-11-13 14:36:48 >>>> generating final output...
FASTA index file /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Adigitifera_genome.fna.fai created.
        >>> Number of toxin loci identified in the genome:
                ACRMI -> 1
                TOTAL -> 1
2025-11-13 14:36:49 >>>> comparing annotated toxins to the uniprot/toxprot database...
        Uniprot/Toxprot database -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta
        >>> Check the uniprot/toxprot report: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Adig/uniprot_report.txt
2025-11-13 14:36:49 >>>> ToxCodAn-Genome finished!
        >>> Check the final annotation output: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Adig/toxin_annotation.gtf

## A tenuis 
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/aten.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Aten -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.



# Expanding database 
## P aus 
### Example Entrez Direct query for cnidarian toxins
esearch -db nucleotide -query "toxin[Title] AND Cnidaria[Organism]" | \
  efetch -format fasta > ncbi_cnidarian_toxins.fasta

#combine all databases together
cat anthozoan_db_gen.fasta ncbi_cnidarian_toxins.fasta > expanded_toxin_db.fasta

sed '/^>/! s/R/A/g' expanded_toxin_db.fasta > expanded_toxin_db_cleaned.fasta


ok lets try again
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/expanded_toxin_db_cleaned.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Paus_expanded -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10

2025-11-13 15:11:18 >>>> searching for toxin regions in genome using the database file...
2025-11-13 15:11:30 >>>> retrieving and filtering matched regions...
2025-11-13 15:11:30 >>>> processing overlapped regions...
2025-11-13 15:11:30 >>>> generating final output...
FASTA index file /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna.fai created.
        >>> Number of toxin loci identified in the genome:
                TOTAL -> 0
2025-11-13 15:11:32 >>>> comparing annotated toxins to the uniprot/toxprot database...
        Uniprot/Toxprot database -> /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta
Warning: [blastx] Query is Empty!
        >>> Check the uniprot/toxprot report: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Paus_expanded/uniprot_report.txt
2025-11-13 15:11:33 >>>> ToxCodAn-Genome finished!
        >>> Check the final annotation output: /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Paus_expanded/toxin_annotation.gtf

  ## P.rus
  python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Prus_genome.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/expanded_toxin_db_cleaned.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results/Prus_expanded -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 20 -l 20 -s 20 -c 10
