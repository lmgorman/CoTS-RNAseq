## Uploading Raw Data to NCBI - 10.03.25
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
Connection terminated after about 20 minutes of trying to move 211R_R1 so only partially moved over in both folders (CoTS_gorman; CoTS_gorman2)

Folder with files (CoTS_gorman; CoTS_gorman2) on NCBI server for 30 days before it will be removed = 09.04.25 - need to upload and submit all data before this date


# 11.03.25
I am going to try and upload all files individually (using just the 'put' command not 'mput') to original 'CoTS_gorman' folder


Navigated to NCBI source folder 
```
cd uploads/kiwielemgee_gmail.com_f7k2V6eS
cd CoTS_gorman
put 218R_R1_001.fastq.gz
```
Completed at 10.51am

```
Transfer complete
5619521752 bytes sent in 08:31 (10.47 MiB/s)
```
Started 218R_R2_001.fastq.gz at 10:52am:
```
put 218R_R2_001.fastq.gz
```
