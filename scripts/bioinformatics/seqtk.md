## Downloading the seqtk software to trim huge raw 211R files (R1 and R2) to subset of sequences
https://github.com/lh3/seqtk 
- seqtk already a module on Unity cluster so no need to download from github

## Make script for seqtk
```
cd /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/scripts
nano seqtk.sh
```
Run script on trimmed 211R files as it saves time vs running on raw data:
```
#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=24:00:00  # Job time limit
#SBATCH -o slurm-seqtk_211R.out  # %j = job ID
#SBATCH -e slurm-seqtk_211R.err  # %j = job ID
#SBATCH -D /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/trimmed_data
```
## Load seqtk module
```
module load seqtk
```

## Subsample 90,000,000 paired reads from the 211R sample
```
seqtk sample -s100 211R_R1_001.fastq 90000000 > sub1_211R_R1.fq
seqtk sample -s100 211R_R2_001.fastq 90000000 > sub2_211R_R2.fq
```

Kept random seed at 100 "-s100" like on documentation

New sample names:
- sub1_211R_R1.fq
- sub2_211R_R2.fq
