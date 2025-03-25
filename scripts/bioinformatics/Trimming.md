Make folders (revise all file paths with real names)
```
pwd
#use this to see what folder you are in and replace directories 

cd /data/gorman/project

mkdir raw-files
mkdir trimmed-files
mkdir scripts
```

Symblic link files from Hollie's folder to your folder 
```
ln -s /project/pi_hputnam_uri_edu/20250107_COTS_LG/*.fastq.gz /data/gorman/project/raw-files

#hard copy
cp /project/pi_hputnam_uri_edu/20250107_COTS_LG/*.fastq.gz /data/gorman/project/raw-files
```

Make a new script document

```
cd /data/gorman/project/scripts/

nano trim-single-file.sh
```

Copy this information into the script document to analyze ONE file. 
```
#!/bin/bash
#SBATCH --job-name=trim-single-file
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=300G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=48:00:00  # Job time limit
#SBATCH -o trim-single-file-%j.out
#SBATCH -e trim-single-file-%j.error
#SBATCH -D /data/gorman/project/raw-files

#load modules
echo "Loading programs" $(date)
module load uri/main
module load fastp/0.23.2-GCC-11.2.0
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b
echo "Starting read trimming." $(date)

#keep this to run on one file

fastp --in1 /data/gorman/project/raw-files/549R_R1_001.fastq.gz  \
--in2 /data/gorman/project/raw-files/549R_R2_001.fastq.gz  \
--out1 /data/gorman/project/trimmed-files/549R_R1_001_trimmed.fastq.gz  \
--out2 /data/gorman/project/trimmed-files/549R_R2_001_trimmed.fastq.gz  \
--detect_adapter_for_pe \
--qualified_quality_phred 20  \
--unqualified_percent_limit 10  \
--cut_right cut_right_window_size 5 cut_right_mean_quality 20 \
-h /data/gorman/project/trimmed-files/549R.fastp.html \
-j /data/gorman/project/trimmed-files/549R.fastp.json

echo "All files trimmed." $(date)
```

Submit the job. 
```
sbatch trim-single-file.sh
```

If you want to run for ALL files, do this: 

```
cd /data/gorman/project/scripts/

nano trim-all-files.sh
```

```
#!/bin/bash
#SBATCH --job-name=trim-all-files
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=300G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=48:00:00  # Job time limit
#SBATCH -o trim-all-files-%j.out
#SBATCH -e trim-all-files-%j.error
#SBATCH -D /data/gorman/project/raw-files

#load modules
echo "Loading programs" $(date)
module load uri/main
module load fastp/0.23.2-GCC-11.2.0
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b
echo "Starting read trimming." $(date)

#keep this to run on all files

array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/gorman/project/trimmed-files/trim.${i} \
        --out2 /data/gorman/project/trimmed-files/trim.$(echo ${i}|sed s/_R1/_R2/) \
		 --detect_adapter_for_pe \
		 --qualified_quality_phred 20  \
		 --unqualified_percent_limit 10  \
		 --cut_right cut_right_window_size 5 cut_right_mean_quality 20 \
		 -h /data/gorman/project/trimmed-files/${i}.fastp.html \
		 -j /data/gorman/project/trimmed-files/${i}.fastp.json		
done

echo "All files trimmed." $(date)
```

Exit nano with control+X  
Then click Y to save   

Submit the job. 
```
sbatch trim-all-files.sh
```

Check the status of the job
```
squeue -u insert-user-name 
```

Cancel all your jobs
```
scancel -u insert-user-name 
```

Cancel one job
```
scancel job-number 
```

You should now see all trim files in your trimmed sequence folder.  
