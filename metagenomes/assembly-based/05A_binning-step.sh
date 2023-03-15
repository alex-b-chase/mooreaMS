#!/bin/bash
#PBS -q home-sio
#PBS -N binMAG_meta
#PBS -l nodes=1:ppn=36:mem768
#PBS -l walltime=38:00:00
#PBS -A jensen-group

module load samtools/1.9
module load hmmer/3.2.1 

source /home/abchase/.bashrc

# run binning on the samples for both maxbin and metabat

# need to generate abundance and depth information on the coverage profiles
# first for metabat

BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly/bowtie
ASSEMBLY=merged.assemblies.s99.fasta

MINLEN=5000

cd $BASEDIR
rm -rf $BASEDIR/metabat

for f in *bowtie2.bam 
do
	newname=${f%.bowtie2.bam}
	samtools sort $f --threads 36 > ${newname}.sorted.bam
done

module load gnu
export CC=/opt/gnu/gcc/bin/gcc
export CXX=/opt/gnu/gcc/bin/g++

jgi_summarize_bam_contig_depths --outputDepth metabat_depth.txt *.sorted.bam

metabat2 -i $ASSEMBLY -a metabat_depth.txt -o metabat/bin \
-m $MINLEN -t 36 -v

echo "metaBAT2 finished successfully, and found $(ls -l metabat | grep .fa | wc -l) bins!"

conda activate py36

# check the binned MAGs with checkm for completeness and purity
cd $BASEDIR
rm -rf metabat/checkm
rm -f metabat.checkm.txt

checkm lineage_wf -f metabat.checkm.txt -t 36 -x .fa metabat/ metabat/checkm

