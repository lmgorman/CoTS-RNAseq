#!/bin/bash
#SBATCH --job-name=truncate_ipr_xml
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=02:00:00
#SBATCH --output=truncate_ipr_xml_%j.log
#SBATCH --error=truncate_ipr_xml_%j.err

# Load Python module if needed
module load python/3.8

# Activate virtualenv if you use one
# source /path/to/venv/bin/activate

# Run your Python script
python /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/truncate_ipr_xml.py \
       /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/iprscan.xml \
       /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/seqid_mapping.tsv \
       /scratch3/workspace/lucy_gorman_uri_edu-lucyscratch/iprscan_truncated.xml
