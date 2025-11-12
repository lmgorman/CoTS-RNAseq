/work/pi_hputnam_uri_edu/lgorman/.conda/envs
active environment : None shell level : 0 user config file : /home/lucy_gorman_uri_edu/.condarc populated config files : /modules/opt/linux-ubuntu24.04-x86_64/miniforge3/24.7.1/.condarc conda version : 24.7.1 conda-build version : not installed python version : 3.12.6.final.0 solver : libmamba (default) virtual packages : __archspec=1=x86_64_v4 __conda=24.7.1=0 __glibc=2.39=0 __linux=6.8.0=0 __unix=0=0 base environment : /modules/opt/linux-ubuntu24.04-x86_64/miniforge3/24.7.1 (read only) conda av data dir : /modules/opt/linux-ubuntu24.04-x86_64/miniforge3/24.7.1/etc/conda conda av metadata url : None channel URLs : https://conda.anaconda.org/conda-forge/linux-64 https://conda.anaconda.org/conda-forge/noarch package cache : /work/pi_hputnam_uri_edu/lgorman/.conda/pkgs envs directories : /home/lucy_gorman_uri_edu/.conda/envs /modules/apps/conda-environments /modules/opt/linux-ubuntu24.04-x86_64/miniforge3/24.7.1/envs platform : linux-64 user-agent : conda/24.7.1 requests/2.32.3 CPython/3.12.6 Linux/6.8.0-62-generic ubuntu/24.04.3 glibc/2.39 solver/libmamba conda-libmamba-so lver/24.7.0 libmambapy/1.5.9 UID:GID : 33364:33364 netrc file : None offline mode : False/home/lucy_gorman_uri_edu/.condarc



#!/bin/bash
#SBATCH --job-name=toxcodan-paus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --constraint=avx512
#SBATCH -p uri-cpu
#SBATCH -o /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/logs/toxcodan-paus-%j.out
#SBATCH -e /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/logs/toxcodan-paus-%j.error

# -----------------------------
# Load modules
# -----------------------------
module load uri/main
module load conda/latest

# Activate environment
conda activate /work/pi_hputnam_uri_edu/lgorman/.conda/envs/ToxcodanGenome


# Run script
python /work/pi_hputnam_uri_edu/lgorman/ToxCodAn-Genome/bin/toxcodan-genome.py \
    -g /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/Paus_genomic.fna \
    -d /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/genome/custom_toxin_database.fasta \
    -o /work/pi_hputnam_uri_edu/20250107_COTS_LG/ToxCodAn/results \
    -c 10
