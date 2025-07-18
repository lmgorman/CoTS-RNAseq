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



