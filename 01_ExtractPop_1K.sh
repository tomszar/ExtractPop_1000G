#!/bin/bash

#This script will download the 1K files, it will keep only SNPs, and only the individuals from EURlist
#Usage:
#    - Place this script in a folder in your work directory
#    - To run it in background and output to a log file run the following code in the console:
#        $ chmod +x 01_ExtractPop_1K.sh
#        $ ./01_ExtractPop_1K.sh > 01_ExtractPop_1K.log 2>&1 &

#Getting the path to this folder
thisdir=$(pwd)

#Download 1000G files if not there and the hg19.fa
cd ~/scratch
for chr in {1..23}
do
    if [ $chr -eq 23 ]; then
        if [ ! -f ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz ]; then
            wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
        fi
        if [ ! -f ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi ]; then
            wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi
        fi
    else
        if [ ! -f ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ]; then
            wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
        fi
        if [ ! -f ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ]; then
            wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
        fi
    fi
done

if [ ! -f hs37d5.fa  ]; then
    wget -q ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    gunzip hs37d5.fa.gz
fi

#Move to this directory
cd $thisdir

for chr in {1..23}
do
    if [ $chr -eq 23 ]; then
        chr=X
    fi
    echo "Writing job for chromosome $chr"
    echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l pmem=32gb
#PBS -A jlt22_b_g_sc_default #jlt22_b_g_sc_default or open
#PBS -j oe

#Moving to working directory
cd ${thisdir}

#Activate conda env
conda activate vcft

#Keep EUR individuals
vcftools --gzvcf ~/scratch/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --remove-indels --keep EURlist.txt --recode --recode-INFO-all --out ~/scratch/EUR_1K_temp_chr${chr}

#Make sure to retain only biallelic SNPs, and remove duplicates
bcftools view -m2 -M2 -v snps ~/scratch/EUR_1K_temp_chr${chr}.recode.vcf | bcftools norm -d none -f ~/scratch/hs37d5.fa | bcftools view -q 0.01:minor | gzip  > ~/scratch/EUR_1K_chr${chr}.vcf.gz

#Not used commands --output-type z --output ~/scratch/EUR_1K_chr${chr}.vcf.gz
#bcftools view EUR_1K_chr22.vcf.gz | awk '/^#/{print}; !/^#/{if (!uniq[\$3]++) print}' | gzip  > Test

rm ~/scratch/EUR_1K_temp_chr${chr}.recode.vcf" >> job_${chr}.pbs

    echo "Submitting job_${chr}.pbs"
    qsub job_${chr}.pbs
    echo "Waiting 1s for next chromosome..."
    sleep 1s 
    echo ""

done
