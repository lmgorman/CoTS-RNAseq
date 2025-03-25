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
Transfer complete
6089363269 bytes sent in 11:08 (8.68 MiB/s)
```

Started 336R_R1_001.fastq.gz at 17:27pm
```
put 336R_R1_001.fastq.gz
Transfer complete
8030925877 bytes sent in 13:31 (9.44 MiB/s)
```

Started 336R_R2_001.fastq.gz at 17:41pm
```
put 336R_R2_001.fastq.gz
Transfer complete
8499368533 bytes sent in 12:57 (10.42 MiB/s)
```

Started 34R_R1_001.fastq.gz at 17:55pm
```
put 34R_R1_001.fastq.gz
Transfer complete
3314376167 bytes sent in 04:49 (10.92 MiB/s)
```
Started 34R_R1_001.fastq.gz at 18:00pm
```
put 34R_R2_001.fastq.gz
Terminal lost connection
```

## 12.03.25
So today I am trying on WindowsPowerShell as I noticed Unity OnDemand had a time out of 15 minutes at top of the On Demand Server - this may have been why my connection was terminated

Started by trying biggest file 211R_R2_001.fastq.gz at 10:47am
```
put 211R_R2_001.fastq.gz
Transfer complete
46647610132 bytes sent in  1:06:30 (11.14 MiB/s)
```

$\color{red}{\textsf{Ok OpenDemand Unity cuts off after 15 minutes so use powershell!}}$ 

Started 34R_R2_001.fastq.gz at 11:54am
```
put 34R_R2_001.fastq.gz
Transfer complete
3481257485 bytes sent in 05:00 (11.04 MiB/s)
```

Started  370R_R1_001.fastq.gz at 12:04pm
```
put  370R_R1_001.fastq.gz
Transfer complete
6449298320 bytes sent in 10:11 (10.05 MiB/s)
```
Started  370R_R2_001.fastq.gz at 12:14pm
```
put  370R_R2_001.fastq.gz
Transfer complete
6806173546 bytes sent in 09:44 (11.10 MiB/s)
```
Started 380R_R1_001.fastq.gz at 12:27pm
```
put  380R_R1_001.fastq.gz
Transfer complete
5935710125 bytes sent in 09:21 (10.07 MiB/s)
```
Started 380R_R2_001.fastq.gz at 12:37pm
```
put  380R_R2_001.fastq.gz
Transfer complete
6182693572 bytes sent in 09:12 (10.66 MiB/s)
```

Started 410R_R1_001.fastq.gz at 12:47pm
```
put 410R_R1_001.fastq.gz
Transfer complete
7446617126 bytes sent in 11:23 (10.39 MiB/s)
```

Started 410R_R2_001.fastq.gz at 12:58pm
```
put 410R_R2_001.fastq.gz
Transfer complete
7700205975 bytes
```

Started 414R_R1_001.fastq.gz at 14:26pm
```
put 414R_R1_001.fastq.gz
Transfer complete
10039996557 bytes sent in 16:22 (9.74 MiB/s)
```
Started 414R_R2_001.fastq.gz at 14:42pm
```
put 414R_R2_001.fastq.gz
Transfer complete
10557980635 bytes sent in 17:28 (9.60 MiB/s)
```
Started 419R_R1_001.fastq.gz at 15:00pm
```
put 419R_R1_001.fastq.gz
Transfer complete
7010002584 bytes sent in 11:35 (9.60 MiB/s)
```

Started 419R_R2_001.fastq.gz at 16:05pm
```
put 419R_R2_001.fastq.gz
Transfer complete
7263856969 bytes sent in 11:35 (9.96 MiB/s)
```
Started 43R_R1_001.fastq.gz at 16:17pm
```
put 43R_R1_001.fastq.gz
Transfer complete
6120591622 bytes sent in 09:04 (10.71 MiB/s)
```

Started 43R_R2_001.fastq.gz at 16:26pm
```
put 43R_R2_001.fastq.gz
Transfer complete
6390333397 bytes sent in 09:59 (10.16 MiB/s)
```
Started 468R_R1_001.fastq.gz at 16:36pm
```
put 468R_R1_001.fastq.gz
Transfer complete
6371225058 bytes sent in 10:07 (10.00 MiB/s)
```
Started 468R_R2_001.fastq.gz at 16:47pm
```
put 468R_R2_001.fastq.gz
Transfer complete
6572885885 bytes sent in 13:25 (7.78 MiB/s)
```
Started 497R_R1_001.fastq.gz at 17:00pm
```
put 497R_R1_001.fastq.gz
Transfer complete
```
Started 497R_R2_001.fastq.gz at 17:12pm
```
put 497R_R2_001.fastq.gz
Transfer complete
6162695969 bytes sent in 09:54 (9.89 MiB/s)
```

Started 512R_R1_001.fastq.gz at 17:22pm
```
put 512R_R1_001.fastq.gz
Transfer complete
7664953698 bytes sent in 12:04 (10.09 MiB/s)
```

## 13.03.25
Started 512R_R2_001.fastq.gz at 10:59am
```
put 512R_R2_001.fastq.gz
Transfer complete
7988983277 bytes sent in 11:58 (10.60 MiB/s)
```
Started 549R_R1_001.fastq.gz at 11:11am
```
put 549R_R1_001.fastq.gz
Transfer complete
5206074821 bytes sent in 07:50 (10.55 MiB/s)
```
Started 549R_R2_001.fastq.gz at 11:19am
```
put 549R_R2_001.fastq.gz
Transfer complete
5453582556 bytes sent in 07:35 (11.40 MiB/s)
```
Started 568R_R1_001.fastq.gz at 11:26am
```
put 568R_R1_001.fastq.gz
Transfer complete
6391960484 bytes sent in 11:00 (9.22 MiB/s)
```
Started 568R_R2_001.fastq.gz at 11:37am
```
put 568R_R2_001.fastq.gz
Transfer complete
6667686438 bytes sent in 10:15 (10.32 MiB/s)
```
Started 571R_R1_001.fastq.gz
```
put 571R_R1_001.fastq.gz
Transfer complete
7136928118 bytes sent in 10:46 (10.52 MiB/s)
```
Started 571R_R2_001.fastq.gz
```
put 571R_R2_001.fastq.gz
Transfer complete
7506695130 bytes sent in 11:53 (10.02 MiB/s)
```
Started 581R_R1_001.fastq.gz at 12:16pm
```
put 581R_R1_001.fastq.gz
Transfer complete
5901914637 bytes sent in 09:13 (10.16 MiB/s)
```

Started 581R_R2_001.fastq.gz at 12:25pm
```
put 581R_R2_001.fastq.gz
Transfer complete
6194699006 bytes sent in 08:49 (11.16 MiB/s)
```
Started 586R_R1_001.fastq.gz at 12:34pm
```
put 586R_R1_001.fastq.gz
Transfer complete
5652557770 bytes sent in 07:47 (11.52 MiB/s)
```
Started 586R_R2_001.fastq.gz at 12:34pm
```
put 586R_R2_001.fastq.gz
Transfer complete
5881746435 bytes sent in 10:44 (8.69 MiB/s)
```
Started 61R_R1_001.fastq.gz at 12:52pm
```
put 61R_R1_001.fastq.gz
Transfer complete
5666923581 bytes sent in 09:01 (9.97 MiB/s)
```
Started 61R_R2_001.fastq.gz at 13:03pm
```
put 61R_R2_001.fastq.gz
Transfer complete
5954980341 bytes sent in 09:39 (9.79 MiB/s)
```
Started 71R_R1_001.fastq.gz at 13:11pm
```
put 71R_R1_001.fastq.gz
Transfer complete
5823011337 bytes sent in 09:21 (9.89 MiB/s)
```
Started 71R_R2_001.fastq.gz at 13:21pm
```
put 71R_R2_001.fastq.gz
Transfer complete
6167294038 bytes sent in 09:45 (10.03 MiB/s)
```
Started 76R_R1_001.fastq.gz at 13:31pm
```
put 76R_R1_001.fastq.gz
Transfer complete
6487379605 bytes sent in 10:03 (10.24 MiB/s)
```
Started 76R_R2_001.fastq.gz at 13:41pm
```
put 76R_R2_001.fastq.gz
Transfer complete
6912449754 bytes sent in 10:19 (10.63 MiB/s)
```
Started 82R_R1_001.fastq.gz at 13:51pm
```
put 82R_R1_001.fastq.gz
Transfer complete
5439924597 bytes sent in 08:22 (10.32 MiB/s)
```
Started 82R_R2_001.fastq.gz at 14:00pm
```
put 82R_R2_001.fastq.gz
Transfer complete
5804378836 bytes sent in 08:01 (11.49 MiB/s)
```
Started 86R_R1_001.fastq.gz at 14:08pm
```
put 86R_R1_001.fastq.gz
Transfer complete
5191500924 bytes sent in 08:02 (10.25 MiB/s)
```
Started 86R_R2_001.fastq.gz at 14:16pm
```
put 86R_R2_001.fastq.gz
Transfer complete
5492476297 bytes sent in 10:00 (8.71 MiB/s)
```

| File name          | Original file size | File size in NCBI | Date uploaded to NCBI |
| ------------------ | ------------------ | ----------------- | --------------------- |
| 16_R1_001.fastq.gz | 6.9G | 6.9G | 10.03.25 |
| 16_R2_001.fastq.gz | 7.2G | 7.2G | 10.03.25 |
| 211R_R1_001.fastq.gz | 42G | 41.2G | 11.03.25 |
| 211R_R2_001.fastq.gz | 44G | | 12.03.25 |
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
| 244R_R2_001.fastq.gz | 6.2G | | 11.03.25 |
| 253R_R1_001.fastq.gz | 6.0G | | 11.03.25 |
| 253R_R2_001.fastq.gz | 6.3G | | 11.03.25 |
| 321RA_R1_001.fastq.gz | 5.6G | | 11.03.25 |
| 321RA_R2_001.fastq.gz | 5.9G | | 11.03.25 |
| 331RA_R1_001.fastq.gz | 5.5G | | 11.03.25 |
| 331RA_R2_001.fastq.gz | 5.7G | | 11.03.25 |
| 336R_R1_001.fastq.gz | 7.5G | | 11.03.25 |
| 336R_R2_001.fastq.gz | 8.0G | | 11.03.25 |
| 34R_R1_001.fastq.gz | 3.1G | | 11.03.25 |
| 34R_R2_001.fastq.gz | 3.3G | | 12.03.25 |
| 370R_R1_001.fastq.gz | 6.1G | | 12.03.25 |
| 370R_R2_001.fastq.gz | 6.4G | | 12.03.25 |
| 380R_R1_001.fastq.gz | 5.6G | | 12.03.25 |
| 380R_R2_001.fastq.gz | 5.8G | | 12.03.25 |
| 410R_R1_001.fastq.gz | 7.0G | | 12.03.25 |
| 410R_R2_001.fastq.gz | 7.2G | | 12.03.25 |
| 414R_R1_001.fastq.gz | 9.4G | | 12.03.25 | 
| 414R_R2_001.fastq.gz | 9.9G | | 12.03.25 | 
| 419R_R1_001.fastq.gz | 6.6G | | 12.03.25 |
| 419R_R2_001.fastq.gz | 6.8G | | 12.03.25 |
| 43R_R1_001.fastq.gz | 5.8G | | 12.03.25 |
| 43R_R2_001.fastq.gz | 6.0G | | 12.03.25 |
| 468R_R1_001.fastq.gz | 6.0G | | 12.03.25 |
| 468R_R2_001.fastq.gz | 6.2G | | 12.03.25 |
| 497R_R1_001.fastq.gz | 5.5G | | 12.03.25 |
| 497R_R2_001.fastq.gz | 5.8G | | 12.03.25 |
| 512R_R1_001.fastq.gz | 7.2G | | 12.03.25 |
| 512R_R2_001.fastq.gz | 7.5G | | 13.03.25 |
| 549R_R1_001.fastq.gz | 4.9G | | 13.03.25 |
| 549R_R2_001.fastq.gz | 5.1G | | 13.03.25 |
| 568R_R1_001.fastq.gz | 6.0G | | 13.03.25 |
| 568R_R2_001.fastq.gz | 6.3G | | 13.03.25 |
| 571R_R1_001.fastq.gz | 6.7G | | 13.03.25 |
| 571R_R2_001.fastq.gz | 7.0G | | 13.03.25 |
| 581R_R1_001.fastq.gz | 5.5G | | 13.03.25 |
| 581R_R2_001.fastq.gz | 5.8G | | 13.03.25 |
| 586R_R1_001.fastq.gz | 5.3G | | 13.03.25 |
| 586R_R2_001.fastq.gz | 5.5G | | 13.03.25 |
| 61R_R1_001.fastq.gz  | 5.3G | | 13.03.25 |
| 61R_R2_001.fastq.gz  | 5.6G | | 13.03.25 |
| 71R_R1_001.fastq.gz  | 5.5G | | 13.03.25 |
| 71R_R2_001.fastq.gz  | 5.8G | | 13.03.25 |
| 76R_R1_001.fastq.gz  | 6.1G | | 13.03.25 |
| 76R_R2_001.fastq.gz  | 6.5G | | 13.03.25 |
| 82R_R1_001.fastq.gz  | 5.1G | | 13.03.25 |
| 82R_R2_001.fastq.gz  | 5.5G | | 13.03.25 |
| 86R_R1_001.fastq.gz  | 4.9G | | 13.03.25 |
| 86R_R2_001.fastq.gz  | 5.2G | | 13.03.25 |

Total GB of folder should be around 456G
```
ls -lh
-rw-rw-r--   1 subftp   submit   7391768242 Mar 10 15:43 16R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7696270619 Mar 10 15:57 16R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   44198100612 Mar 11 13:14 211R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   46647610132 Mar 12 11:53 211R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5619521752 Mar 11 10:51 218R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6013697096 Mar 11 11:01 218R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   2948199684 Mar 11 11:06 225R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   3060671251 Mar 11 11:11 225R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7235529321 Mar 11 11:23 227R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7513266707 Mar 11 11:35 227R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6397379460 Mar 11 11:46 235R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6832723749 Mar 11 11:58 235R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6998357089 Mar 11 12:08 236R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7347856690 Mar 11 13:25 236R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6228927206 Mar 11 13:34 244R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6603779223 Mar 11 16:22 244R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6368047803 Mar 11 16:31 253R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6723438891 Mar 11 16:42 253R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6001367402 Mar 11 16:52 321RA_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6269993774 Mar 11 17:02 321RA_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5813803998 Mar 11 17:16 331RA_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6089363269 Mar 11 17:27 331RA_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   8030925877 Mar 11 17:41 336R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   8499368533 Mar 11 17:54 336R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   3314376167 Mar 11 18:00 34R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   3481257485 Mar 12 11:58 34R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6449298320 Mar 12 12:13 370R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6806173546 Mar 12 12:23 370R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5935710125 Mar 12 12:36 380R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6182693572 Mar 12 12:46 380R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7446617126 Mar 12 12:58 410R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7700205975 Mar 12 13:09 410R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   10039996557 Mar 12 14:42 414R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   10557980635 Mar 12 14:59 414R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7010002584 Mar 12 15:11 419R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7263856969 Mar 12 16:17 419R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6120591622 Mar 12 16:26 43R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6390333397 Mar 12 16:36 43R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6371225058 Mar 12 16:46 468R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6572885885 Mar 12 17:00 468R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5867836245 Mar 12 17:09 497R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6162695969 Mar 12 17:22 497R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7664953698 Mar 12 17:34 512R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7988983277 Mar 13 11:10 512R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5206074821 Mar 13 11:18 549R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5453582556 Mar 13 11:26 549R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6391960484 Mar 13 11:37 568R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6667686438 Mar 13 11:48 568R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7136928118 Mar 13 12:03 571R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   7506695130 Mar 13 12:15 571R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5901914637 Mar 13 12:25 581R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6194699006 Mar 13 12:33 581R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5652557770 Mar 13 12:41 586R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5881746435 Mar 13 12:52 586R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5666923581 Mar 13 13:01 61R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5954980341 Mar 13 13:11 61R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5823011337 Mar 13 13:21 71R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6167294038 Mar 13 13:30 71R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6487379605 Mar 13 13:41 76R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   6912449754 Mar 13 13:51 76R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5439924597 Mar 13 14:00 82R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5804378836 Mar 13 14:08 82R_R2_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5191500924 Mar 13 14:16 86R_R1_001.fastq.gz
-rw-rw-r--   1 subftp   submit   5492476297 Mar 13 14:26 86R_R2_001.fastq.gz
```
Total files = 64 (32 samples x 2 reads _per_ sample)
Total file size was 455.2GB
Submitted to NCBI for processing on 13.03.25 at 14:44pm
