#!/bin/bash

#This script will download the 1K files, it will keep only SNPs, and only the individuals from EURlist
#Usage:
#    - Place this script in a folder in your work directory
#    - To run it in background and output to a log file run the following code in the console:
#        $ chmod +x 01_ExtractPop_1K.sh
#        $ ./01_ExtractPop_1K.sh > 01_ExtractPop_1K.log 2>&1 &

#Getting the path to this folder
thisdir=$(pwd)

#Download 1000G files if not there
cd ~/scratch
for chr in {1..22}
do
    if [ ! -f ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ]; then
        wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    fi
    if [ ! -f ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ]; then
        wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
    fi
done

#Move to this directory
cd $thisdir

for chr in {1..22}
do
    echo "Writing job for $chr"
    echo "#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l pmem=16gb
#PBS -A jlt22_b_g_sc_default #jlt22_b_g_sc_default or open
#PBS -j oe

#Moving to working directory
cd ${thisdir}

#Activate conda env
conda activate vcft

vcftools --gzvcf ~/scratch/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --remove-indels --keep EURlist.txt --out ~/scratch/EUR_1K_chr${chr} --plink" >> job_${chr}.pbs

    echo "Submitting job_${chr}.pbs"
    qsub job_${chr}.pbs
    echo "Waiting 5s for next chromosome..."
    sleep 5s 

done
