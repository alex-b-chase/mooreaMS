#!/bin/bash
#PBS -q home-sio
#PBS -N binMAG_maxbin
#PBS -l nodes=1:ppn=16
#PBS -l walltime=18:00:00
#PBS -A jensen-group

module load samtools/1.9
module load hmmer
module load idba-ud/1.1.1
module load bowtie2/2.3.4.3 

source /home/abchase/.bashrc

# now for maxbin

BASEDIR=/oasis/tscc/scratch/abchase/moorea
ASSEMBLY=merged.assemblies.s99.fasta

MINLEN=2000

cd $BASEDIR

cut -f1,3 metabat_depth.txt | tail -n+2 > maxbin.abund_list.txt

rm -rf $BASEDIR/maxbin
mkdir -p $BASEDIR/maxbin

run_MaxBin.pl -contig $ASSEMBLY -out maxbin/bin -abund maxbin.abund_list.txt \
-min_contig_length $MINLEN -thread 16 -verbose -markerset 40

rm -rf $BASEDIR/maxbin2
mkdir -p $BASEDIR/maxbin2

N=0
for i in $(ls $BASEDIR/maxbin/ | grep bin | grep .fasta)
do
	cp $BASEDIR/maxbin/$i $BASEDIR/maxbin2/bin.${N}.fa
	N=$((N + 1))
done

echo "MaxBin2 finished successfully, and found $(ls -l $BASEDIR/maxbin2 | grep .fa | wc -l) bins!"

conda activate py36

# check the binned MAGs with checkm for completeness and purity
cd $BASEDIR
rm -rf maxbin2/checkm
rm -f metabat.checkm.txt

checkm lineage_wf -f metabat.checkm.txt -t 16 -x .fa maxbin2/ maxbin2/checkm

