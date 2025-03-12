## Coding - Lucy's Intro

Ariana’s terminal command notes 
https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/2023-UW-Software-Carpentry-Workshop-in-Terminal-Git-and-R/ 

## Useful commands
```
ls = lists everything in current directory you are in
ls -h -l = lists everything PLUS size of files
Everything before $ sign tells you where you are
pwd = prints working directory 
cd = change directory “cd mcap-2023”
cd .. = goes back one level in folders
rm = removes a file
mkdir = makes a folder in the current directory
nano = good text editor program
.sh = shell script
module --show_hidden spider fastp = shows whether fastp has been downloaded as a program in the HPC - helps you see whether a program is present in the HPC 

To run a shell script “sbatch align.sh” = submits job you have coded in nano as a shell script 
Always a header at the start of the scripts in the nano text software as the computer reads this to understand how much memory and cpu to dedicate to the jobs
scancel #job = deletes job 
squeue -u username
-D = where to look for all data
-o = output file
-e = error file
less = view file
echo = repeat
exit = exits out of interactive mode
module avail = all programs available on supercomputer
array = make list 
cp = copy
* = Asterix refers to everything 
interactive = mini session with supercomputer
Run test script on one file in interactive mode = lets you check script works before doing all files
```

## Supercomputer header for each shell script
```
“#SBATCH” = language that supercomputer needs to read
e.g. each nano shell script needs this at the beginning:
#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-fastqc_raw.out  # %j = job ID
#SBATCH -e slurm-fastqc_raw.err  # %j = job ID
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/raw_data
```

## Easy example nano shell script
```
#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-fastqc_raw.out  # %j = job ID
#SBATCH -e slurm-fastqc_raw.err  # %j = job ID
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/raw_data

#load modules 
module load uri/main
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b

#run fastqc on raw data
fastqc /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/raw_data/*.fastq.gz -o /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/raw_multiqc/

#generate multiqc report
multiqc /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/raw_multiqc/ --filename multiqc_report_raw.html 

echo "Initial QC of raw seq data complete." $(date)
Then close the shell script and save as "raw_qc.sh"

#Run the job, type:
sbatch raw_qc.sh
```

## File Transfer Protocols
Use the command 
```
ftp
```
and then enter the server to transfer the files to
You can use the command 
```
mput
```
to move files to the server

## Github text editing
https://gist.github.com/luigiMinardi/4574708d404cdf4fe0da7ac6fe2314db#colors = changing colour texts in markdowns

```diff
- Red
+ Green
! Orange
# Grey
@@ Purple @@
```

