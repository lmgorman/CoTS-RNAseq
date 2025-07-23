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

Based on this
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


