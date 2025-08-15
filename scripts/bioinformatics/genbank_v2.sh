#Convert to gff3
gffread /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff -o GCA_964291705.1_jaAcrHyac4.1_genomic_2.gff
#convert to gtf
gffread /work/pi_hputnam_uri_edu/ashuffmyer/cots-gorman/acr_hya_jaAcrHyac4.1/GCA_964291705.1_jaAcrHyac4.1_genomic.gbff -T -o GCA_964291705.1_jaAcrHyac4.1_genomic_2.gtf
