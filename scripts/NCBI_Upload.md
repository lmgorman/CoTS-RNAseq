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

I am uploading the files in sequential order when you list the files in the folder:
```
/project/pi_hputnam_uri_edu/20250107_COTS_LG ls
```
The ONLY exception = I am leaving the 211R_R1 and 211R_R2 until last due to their size (takes around one hour to transfer one file to NCBI)

Navigated to NCBI source folder and started transferring 218R_R1_001.fastq.gz at 10:43am
```
cd uploads/kiwielemgee_gmail.com_f7k2V6eS
cd CoTS_gorman
put 218R_R1_001.fastq.gz
Transfer complete
5619521752 bytes sent in 08:31 (10.47 MiB/s)
```
Started 218R_R2_001.fastq.gz at 10:52am:
```
put 218R_R2_001.fastq.gz
Transfer complete
6013697096 bytes sent in 09:16 (10.30 MiB/s)
```
Started 225R_R1_001.fastq.gz at 11:02am

```
put 225R_R1_001.fastq.gz
Transfer complete
2948199684 bytes sent in 04:32 (10.32 MiB/s)
```

Started 225R_R2_001.fastq.gz at 11:06am
```
put 225R_R2_001.fastq.gz
Transfer complete
3060671251 bytes sent in 04:30 (10.78 MiB/s)
```

Started 227R_R1_001.fastq.gz at 11:12am
```
put 227R_R1_001.fastq.gz
Transfer complete
7235529321 bytes sent in 11:12 (10.25 MiB/s)
```
Started 227R_R2_001.fastq.gz at 11:23am
```
put 227R_R2_001.fastq.gz
Transfer complete
7513266707 bytes sent in 12:14 (9.75 MiB/s)
```
Started 235R_R1_001.fastq.gz at 11:36am
```
put 235R_R1_001.fastq.gz
```
Connection to remote server was terminated during = seems that I can access for one hour

Relogged in:
```
cd /project/pi_hputnam_uri_edu/20250107_COTS_LG
Connected using ftp credentials
cd uploads/kiwielemgee_gmail.com_f7k2V6eS
cd CoTS_gorman
put 235R_R1_001.fastq.gz
Transfer complete
6397379460 bytes sent in 02:24 (42.17 MiB/s)
```
Taking way less time ~only 2 mins? Has this transfer worked correctly? Need to check size of file in folder on NCBI

Started 235R_R2_001.fastq.gz at 11:48am
```
put 235R_R2_001.fastq.gz
Transfer complete
6832723749 bytes sent in 10:51 (10.00 MiB/s)
```
Started 236R_R1_001.fastq.gz at 11:58am
```
put 236R_R1_001.fastq.gz
Transfer complete
6998357089 bytes sent in 09:39 (11.50 MiB/s)
```
Started 211R_R1_001.fastq.gz at 12:11am
```
put 211R_R1_001.fastq.gz
```
Re-uploading the whole file versus the 235R_R1_001.fastq.gz this morning - maybe if you log back in on the same day it carries on where the file left off uploading? Will try this if connection gets terminated during upload....
Connection terminated after 25 minutes; checked and if i log in again straight after it looks to be transferring reads from where it left off (similar to 235R_R1_001.fastq.gz this morning):
only 16 minutes remaining - need to check true file size after upload completed
```
put 211R_R1_001.fastq.gz
Transfer complete
44198100612 bytes sent in 28:50 (24.35 MiB/s)
```
Started 236R_R2_001.fastq.gz at 13:14pm 
```
put 236R_R2_001.fastq.gz
Transfer complete
7347856690 bytes sent in 10:51 
```
Started 244R_R1_001.fastq.gz at 13:25pm
```
put 244R_R1_001.fastq.gz
Transfer complete
6228927206 bytes sent in 09:11 (10.76 MiB/s)
```
Started 211R_R1_001.fastq.gz at 14:55pm
```
put 211R_R2_001.fastq.gz 
```
- did 50/55 minutes and got to ~90% and connection terminated

- Started re running at 15:55pm - can't tell whether its started from the beginning again or just continuing from where it left off....
```
put 211R_R2_001.fastq.gz 
```
I exited because it was just going to do the same again - try and transfer for 50 minutes and then loose conneciton...
*Not completed 211R_R2_001.fastq.gz*

Started 244R_R2_001.fastq.gz at 16:12pm
```
put 244R_R2_001.fastq.gz
Transfer complete
6603779223 bytes sent in 09:39 (10.87 MiB/s)
```
Started 253R_R1_001.fastq.gz at 16:22pm
```
put 253R_R1_001.fastq.gz
Transfer complete
6368047803 bytes sent in 09:47 (10.33 MiB/s)
```
Started 253R_R2_001.fastq.gz at 16:34pm
```
put 253R_R2_001.fastq.gz
Transfer complete
6723438891 bytes sent in 10:30 (10.16 MiB/s)
```
Started 321RA_R1_001.fastq.gz at 16:42pm
```
put 321RA_R1_001.fastq.gz
Transfer complete
6001367402 bytes sent in 10:23 (9.18 MiB/s)
```
Started 321RA_R2_001.fastq.gz at 16:53pm
```
put 321RA_R2_001.fastq.gz
Transfer complete
6269993774 bytes sent in 09:35 (10.39 MiB/s)
```
Started 331RA_R1_001.fastq.gz at 17:06pm
```
put 331RA_R1_001.fastq.gz
Transfer complete
5813803998 bytes sent in 10:08 (9.10 MiB/s)
```
Started 331RA_R2_001.fastq.gz at 17:16pm
```
put 331RA_R2_001.fastq.gz
```


| File name          | Original file size | File size in NCBI | Date uploaded to NCBI |
| ------------------ | ------------------ | ----------------- | --------------------- |
| 16_R1_001.fastq.gz | 6.9G | 6.9G | 10.03.25 |
| 16_R2_001.fastq.gz | 7.2G | 7.2G | 10.03.25 |
| 211R_R1_001.fastq.gz | 42G | 41.2G | 11.03.25 |
| 211R_R2_001.fastq.gz | 44G | | |
| 218R_R1_001.fastq.gz | 5.3G | 5.2G | 11.03.25 |
| 218R_R2_001.fastq.gz | 5.7G | 5.6G | 11.03.25 |
| 225R_R1_001.fastq.gz | 2.8G | | 11.03.25 |
| 225R_R2_001.fastq.gz | 2.9G | | 11.03.25 |
| 227R_R1_001.fastq.gz | 6.8G | | 11.03.25 |
| 227R_R2_001.fastq.gz | 7.0G | | 11.03.25 |
| 235R_R1_001.fastq.gz | 6.0G | | 11.03.25 |
| 235R_R2_001.fastq.gz | 6.4G | | 11.03.25 |
| 236R_R1_001.fastq.gz | 6.6G | | 11.03.25 |
| 236R_R2_001.fastq.gz | 6.9G | | 11.03.25 |
| 244R_R1_001.fastq.gz | 5.9G | | 11.03.25 |
| 244R_R2_001.fastq.gz | | | 11.03.25 |
| 253R_R1_001.fastq.gz | | | 11.03.25 |
| 253R_R2_001.fastq.gz | | | 11.03.25 |
| 321RA_R1_001.fastq.gz | | | 11.03.25 |
| 321RA_R2_001.fastq.gz | | | 11.03.25 |
| 331RA_R1_001.fastq.gz | | | 11.03.25 |
| 331RA_R2_001.fastq.gz | | |
| 336R_R1_001.fastq.gz | | |
