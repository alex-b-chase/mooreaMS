#!/bin/bash
#PBS -q home-sio
#PBS -N collateMAG_dRep
#PBS -l nodes=1:ppn=8:mem384
#PBS -l walltime=17:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4


BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly
OUTDIR=$BASEDIR/bowtie

conda activate checkm

export PATH="$PATH:/home/abchase/software/mash-Linux64-v2.3"
export PATH="$PATH:/home/abchase/software/mummer-4.0.0rc1"
export PATH="$PATH:/home/abchase/software/FastANI-master"
export PATH="$PATH:/home/abchase/software/ANIcalculator_v1"

cd $OUTDIR/combinedMAGs

rm -rf nt_files
rm -rf aa_files
rm -rf origMAG
mkdir -p nt_files
mkdir -p aa_files
mkdir -p origMAG

for f in *.fa
do
	newname=${f%.fa}
	prodigal \
	-i $f -d ${newname}.fna \
	-a ${newname}.faa -m -p meta 

	mv $f origMAG
	mv ${newname}.fna nt_files
	mv ${newname}.faa aa_files
done

rm -rf checkM_out
rm -rf testDrep

checkm lineage_wf \
aa_files checkM_out \
-f checkM_out/metabat.checkm.tsv \
-g -x faa --tab_table \
-t 8 --pplacer_threads 8 

checkm qa \
checkM_out/lineage.ms checkM_out \
-f checkM_out/Chdb.tsv \
-t 8 --tab_table -o 2

echo "genome,completeness,contamination" > genomeInformation.csv
cut -f1,6,7 checkM_out/Chdb.tsv | awk '{$1=$1 ".fa"}1' | sed 's/ /,/g' | grep -v "Completeness" >> genomeInformation.csv

dRep dereplicate \
testDrep -p 8 \
-g origMAG/*.fa --genomeInfo genomeInformation.csv

conda deactivate

cd $OUTDIR/combinedMAGs/testDrep
rm -rf gtdbclass

conda activate gtdbtk

gtdbtk classify_wf \
--genome_dir dereplicated_genomes \
--out_dir gtdbclass \
-x fa --cpus 8 --force 

conda deactivate





