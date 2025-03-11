## Uploading Raw Data to NCBI
Follow steps detailed at: https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md

1. Add biosamples and attributes
[Gorman_CoTS-coral_250225_MIMS.me.host-associated.6.0.xlsx](https://github.com/user-attachments/files/19166056/Gorman_CoTS-coral_250225_MIMS.me.host-associated.6.0.xlsx)
2. Add bioproject for biosamples
3. Add SRA metadata
[Gorman_CoTScoral_SRA_metadata_acc (2).xlsx](https://github.com/user-attachments/files/19166044/Gorman_CoTScoral_SRA_metadata_acc.2.xlsx)
4. Upload data

Navigate to folder with raw data on Unity
```
/project/pi_hputnam_uri_edu/20250107_COTS_LG
```
Establish remote connection using ftp command
```
ftp ftp-private.wip.ncbi.nlm.nih.gov
```
Login using credentials NCBI give you
Now you are in NCBI server
Create folder in NCBI server
```
mkdir CoTS_gorman
```
Navigate to folder in NCBI


Transfer multiple files in your home directory folder (/project/pi_hputnam_uri_edu/20250107_COTS_LG) using 'mput' command 
```
mput *fastq.gz
```
Then the ftp will ask '[anpqy?]'
And you want to type 'a' for all files

This will transfer all the files ending in "fastq.gz" to your NCBI folder

As of 10.03.25
2 samples fully uploaded in NCBI files:
CoTS_gorman (16_R1; 16_R2)
CoTS_gorman2 (16_R1; 16_R2)
