#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/capra/projects/methylation_clocks/Aim_1/MAGENTA_genetics/get_individual_ancestry_loci.log
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00

# This script gathers, for all genotyped loci, the local ancestry of a given individual for each haplotype, as well as the genotype for the individual.

module load CBI
module load bcftools

# Use bcftools query. Get the sample ID, the "ID" of the genotyped locus, the two ancestries, and the genotype.
bcftools query -f '[%SAMPLE\t %ID\t %AN1\t %AN2\t %GT\n]' AA_merged.vcf.gz > AA_individual_ancestry_loci.txt
bcftools query -f '[%SAMPLE\t %ID\t %AN1\t %AN2\t %GT\n]' NHW_merged.vcf.gz > NHW_individual_ancestry_loci.txt
bcftools query -f '[%SAMPLE\t %ID\t %AN1\t %AN2\t %GT\n]' HISPANIC_merged.vcf.gz > HISPANIC_individual_ancestry_loci.txt