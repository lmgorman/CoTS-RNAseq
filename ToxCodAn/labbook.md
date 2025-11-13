# 12.11.25
Created new directory under
`/work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn`


Created new genome directory to save the genomes into
`/work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome`

To use the method from
(https://github.com/pedronachtigall/ToxCodAn-Genome)



Downloaded their anthozoan CDS
https://github.com/pedronachtigall/ToxCodAn-Genome/blob/main/Databases/Anthozoa_db.fasta


# 13.11.25
## P. aus
Ran P.aus genome on default settings
`python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -c 10`

Didnt return any hits


Changed % identity to 40% and minimu  length from 200 to 20
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 40 -l 20 -c 10
ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

Changed gene size:
    python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 40 -l 20 -s 20 -c 10
  ToxCodAn-Genome stopped to run, due to the following issue:
  No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

  ## P.rus
  python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Prus_genome.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 40 -l 20 -s 20 -c 10

  ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.

## P. compressa
  Porites_compressa.fasta
   python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Porites_compressa.fasta -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/anthozoan_db_gen.fasta -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results -u /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta -p 40 -l 20 -s 20 -c 10

   ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.
    
