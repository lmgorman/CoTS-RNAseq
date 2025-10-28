## 22.04.25 - Ran InterPro Scan on Acropora hyacinthus genome
Made script - [interpro2.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/275c71ed56a750e267e850c7af951438dafb736a/scripts/bioinformatics/interpro2.sh) = uses new InterPro scan version
```
Job ID: 32969186
```

## 25.04.25 - Ran InterPro Scan on Porites evermanni genome
Made script - [interpro_pever.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/537510c9e60907224b9f36664a9a88307ce7f927/scripts/bioinformatics/interpro_pever.sh) Started running on 25.04.25
```
Job ID:33156317
squeue -j 33156317
```

Had to clean up protein fasta as it contained “*” - made script called [remove_asterix.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/537510c9e60907224b9f36664a9a88307ce7f927/scripts/bioinformatics/remove_asterix.sh) to do this

Cleaned fasta file now found: 
```
/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/refs/por/Porites_evermanni_v1_clean.annot.pep.fa
```
Done - checked on 26.04.25


## 26.04.25 - Ran eggnog mapper on Acropora hyacinthus & Porites evermanni genomes
Scripts can be found here: 

[scripts/bioinformatics/eggnog_Ahya.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/537510c9e60907224b9f36664a9a88307ce7f927/scripts/bioinformatics/eggnog_Ahya)

[scripts/bioinformatics/eggnog_Pever.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/537510c9e60907224b9f36664a9a88307ce7f927/scripts/bioinformatics/eggnog_Pever.sh)
```
eggnog_Ahya.sh JOB ID: 33386208
eggnog_Pever.sh JOB ID: 33399582
```

## 16.07.25 - Run funannotate predict on Acropora hyacinthus genome
Added —force to script as it was terminating after bad contigs


Reran
```
JOB ID:40032166
squeue -j 40032166
```
Started at 15:15pm


Script can be found here:
[scripts/bioinformatics/funn_pred_ahya.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/04aabf2b71618d285860fcc2b134ea8d8f964345/scripts/bioinformatics/funn_pred_ahya.sh) 

Error file says 20 bad contigs - I just force passed the contigs (not reccommended) probably because I didn’t clean the genome using funannotate clean

## 17.07.25
```
sacct -j 40032166 --format=JobID,State,ExitCode


JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         40032166   uri-cpu funannot lucy_gor  R   21:30:50      1 uri-cpu009
```

The time constraint on the job is 48 hours… change number of days to 5 to shell scripts

**DONE**

To determine quality:
  - Look at number of genes annotated in total
  - How many GO terms and general terms assigned to each gene 
  - What % in total annotated genes
  - more than 50% annotation = good
  - Acropora is well annotated so higher percentage for Acropora vs Porites 
  - Can re run with optimised augustus and compare % of overall genes annotated


Error file said:
```
"Funannotate predict is finished, output files are in the /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/Ahya_funann//predict_results folder"
```
Outfile:
```
2742 completed 0 remaining 
You have 43,700 genes predicted overall.
"mRNA": 38285
"tRNA": 5415
CDS_complete: 37,751
total_exons: 226,958 exons
```

## 18.07.25 - Run funannotate predict on Acropora hyacinthus genome with optimised Augustus
Running optimised prediction (funn_pred_opti_Ahya.sh) using ``--optimise Augustus`` parameter and will compare that to the non-optimised I did yesterday


Script found here:

[scripts/bioinformatics/funn_pred_opti_Ahya.sh ](https://github.com/lmgorman/CoTS-RNAseq/blob/09e8f1ef8ca2a538dfca4293b2b3bbdb6acf3dc7/scripts/bioinformatics/funn_pred_opti_Ahya.sh)

Started funn_pred_opti_Ahya.sh 18.07.25 12:51pm
squeue -j 40082840


Saving the data to my scratch directory so I don't take up more space in the /work/pi_hputnam_uri_edu


Next:
- Decide which prediction you will use optimised augustus or not
- Then proceed to annotation For A.hyacinthus!
- Repeat for P. evermanni
- Will need to clean the P. evermanni genome also with funnannotate clean

## 23.07.25
My structural annotations do not fit the functional annotation as well as Lopez-Nandam - my structural had 0 GO terms and 0 Eggnog terms


Therefore using Lopez Nandam structural annotations and the functional annotation based on this can be found at:
/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/Ahya_ann/functional_annotate_results

Based on this:

**Table 1.** Functional annotation results of *Acropora hyacinthus* genome using funannotate (version 1.8.17).
| Tool/DB |Hits| % of total transcripts (n = 27,110) |
|---------|----|-------------------------------------|
|GO Terms|13,557|~50%|
|InterProScan|18,149|~67%|
|eggNOG|20,372|~75%|
|Pfam|12,804|~47%|
|CAZyme|216|<1% (expected: low)|
|MEROPS|782|~2.9%|
|BUSCO|780|(possibly number of BUSCO genes found)|
|Secretion|0| Didn't run secretion software|

Downloaded this file from unity to local PC using instructions on:
https://docs.unity.uri.edu/documentation/managing-files/filezilla/#set-up-default-remote-directory-in-filezilla 

## Downloaded new Acropora hyacinthus genome
More updated genome with better BUSCO scores - see https://gigabytejournal.com/articles/153 

**Table 2.** Comparison of genome quality (BUSCO scores) between *Acropora hyacinthus* genomes from Lopez-Nandam and Aquatic symbiosis project.
| Database/Study | BUSCO complete | BUSCO fragmented | missing |
|----------------|----------------|------------------|---------|
|Lopez-Nandam| 71.3%| 13.6%| 15.15%|
|Aquatic symbiosis| 96.7% (single copy 95.5%; dual copy 1.2%)| 1.2% |2.2%|

Files available here - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/291/705/GCA_964291705.1_jaAcrHyac4.1/

Uploaded the gbff and fasta to:
```
/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/
```

Might need to re-map all RNASeq data using STAR
First need to convert to gbff to gff3 format


Asked HPC mass to download:
- genome tools to do this
git clone https://github.com/genometools/genometools.git
cd genometools
make
sudo cp bin/gt /usr/local/bin/


and also download 
- genbank_to command in python
https://pypi.org/project/genbank-to/


I also need to convert the genome fasta to protein sequences for InterProScan and eggnog mapper (the genome doesn't have this data)
looks like biopython is a commonly used tool to do this
`from Bio import SeqIO`


Current scripts that need updating based on the protein sequences and the gff3 read format:
- [build_star_index_jaAcrHyac4.1.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/3d98e29810ee6c4c8b61581b6480c12db3780f26/scripts/bioinformatics/build_star_index_jaAcrHyac4.1.sh) = trying to index the new genome and then it can be used for mapping - waiting on converting genbank to gff3
- [map_549R_newAhyagen.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/6ab62eeacc69884b5b311b4182287f228bc7f8f2/scripts/bioinformatics/map_549R_newAhyagen.sh) = testing out mapping my RNASeq 549R data to new A. hyacinthus genome - waiting on script above being successful
- [ahya_NEW_interpro.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/152c4204f68dbb233f5f97e463eb64d76fbe195c/scripts/bioinformatics/ahya_NEW_interpro.sh) = InterProscan data for new A. hyacinthus - waiting on using program to convert DNA to proteins to input into script

## 14.08.25 - Made script to convert genbank to gff3 using biopearl
#Cloned the bioperl git hub repo
```
 cd /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman
 git clone https://github.com/bioperl/bioperl-live.git
```

Ran the following script in /home/lucy_gorman_uri_edu/scripts:
[genbank_to_gff3.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/d78867a608c706ae76f53f76d0b8f50987e57499/scripts/bioinformatics/genbank_to_gff3.sh)
```
squeue -j 41061834
GFF3 saved to /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff.gff
```
Ok now I need to convert the GFF3 to GTF to run in STAR index so:
[build_star_index_jaAcrHyac4.1.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/3d98e29810ee6c4c8b61581b6480c12db3780f26/scripts/bioinformatics/build_star_index_jaAcrHyac4.1.sh)

squeue -j 41062013

Error file in /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/build_star_index_jaAcrHyac4.1-41062013.error:
```
Loading uri version main
NOTE: The modules under this branch will not run on the login node. Use
--constraint=avx512 for sbatch or srun sessions.
Loading gffread version 0.12.7
!!!!! WARNING: Could not move Log.out file from ./Log.out into /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/star_index_jaAcrHyac4.1/Log.out. Will keep ./Log.out


Fatal INPUT FILE error, no exon lines in the GTF file: GCA_964291705.1_jaAcrHyac4.1_genomic.gtf
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

Aug 14 14:00:46 ...... FATAL ERROR, exiting
cp: cannot stat '/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/star_index_jaAcrHyac4.1/*': No such file or directory
build_star_index_jaAcrHyac4.1-41062013.error (END)Loading uri version main
NOTE: The modules under this branch will not run on the login node. Use
--constraint=avx512 for sbatch or srun sessions.
Loading gffread version 0.12.7
!!!!! WARNING: Could not move Log.out file from ./Log.out into /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/star_index_jaAcrHyac4.1/Log.out. Will keep ./Log.out


Fatal INPUT FILE error, no exon lines in the GTF file: GCA_964291705.1_jaAcrHyac4.1_genomic.gtf
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

Aug 14 14:00:46 ...... FATAL ERROR, exiting
cp: cannot stat '/scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/star_index_jaAcrHyac4.1/*': No such file or directory
```

GTF file issue:

The file GCA_964291705.1_jaAcrHyac4.1_genomic.gtf does not contain any lines where the 3rd column is exon.

STAR requires exon features in the GTF file to build the splice junction database

Chat to Ariana tomorrow about editing the file to change this

# 22.10.25
Sent myself the P. evermanni genbank file combined with genome file made in Geneious Prime as a `.gff` file

Upload this to unity and retry the funannotate command on Unity for this genome

The new A. hyacinthus genome doesn't have any annotated mRNA or exon features and you need transcriptomics data to do this - so I cannot index it using STAR

 Ariana said I can try this without structural annotations and just compare % mapped reads to my RNA data 

# 26.10.25
Downloaded the gff file from Geneious of Porites evermanni gff


I FileZilla transferred it to:
/work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/por-ever/

Now will do this and retry funannotate script - funann_pever.sh


file called: Porites_evermanni_v1_CORRECT.gff

Still not great so using agat in bioconda
script = [agat_fix_Pever_gff.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/73908b0ecb13bf3d030bda7ca31a761eefa2d716/scripts/bioinformatics/agat_fix_Pever_gff.sh)

Submitted job - squeue -j 48107324
Worked woo!

Now running funannotate on P. evermanni genome:
script: [funnann_pever.sh](https://github.com/lmgorman/CoTS-RNAseq/blob/9ecbe3c1fee78f156e807449cc1203366653a6b6/scripts/bioinformatics/funann_pever.sh)
squeue -j 48107358

Locus of genbank file crashed funannotate
running new job squeue -j 48153862


Running
awk 'BEGIN{OFS="\t"} $3=="CDS" && $8=="."{$8=0}1' Porites_evermanni_v1_FIXED.short.gff > Porites_evermanni_v1.short.fixphase.gff 

to remove "." in CDS phase
squeue -j 48201829
