#Create conda environment in scratch directory in an interactive session
salloc -p cpu -c 2 --mem=5G --time 02:00:00
conda activate ToxcodanGenome

Create the environment:
conda create -n ToxcodanGenome -c bioconda python biopython pandas blast exonerate miniprot gffread hisat2 samtools stringtie trinity spades
Git clone the ToxCodAn-Genome repository and add the bin to your PATH:
git clone https://github.com/pedronachtigall/ToxCodAn-Genome.git
echo "export PATH=$PATH:$PWD/ToxCodAn-Genome/bin/" >> ~/.bashrc
Replace ~/.bashrc to ~/.bash_profile if needed.
It may be needed to apply "execution permission" to all bin executables in "ToxCodAn-Genome/bin/":
chmod +x ToxCodAn-Genome/bin/*
Then, run ToxCodAn-Genome as described in the "Usage" section.
To activate the environment to run ToxCodAn-Genome just use the command: conda activate ToxcodanGenome
To deactivate the environment just use the command: conda deactivate
