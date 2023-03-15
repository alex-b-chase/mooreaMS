#!/bin/bash

BASEDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/binned_MAGs/dRep
GENOMEDIR=$BASEDIR/dereplicated_genomes
OUTDIR=$BASEDIR/binMAGanalysis

rm -rf $OUTDIR
mkdir -p $OUTDIR

cp $GENOMEDIR/*.fa $OUTDIR/
cd $OUTDIR

rm -rf $OUTDIR/prokka_annotation
rm -rf $OUTDIR/antismash
mkdir -p $OUTDIR/prokka_annotation
mkdir -p $OUTDIR/antismash


for genome in *.fa
do
	cd $OUTDIR
	genomename=${genome%.fa}

	cp $genome $OUTDIR/prokka_annotation/${genomename}.fna

	cd $OUTDIR/prokka_annotation
	prokka --outdir ${genomename} \
	--prefix ${genomename} --kingdom Bacteria \
	--notrna --quiet \
	${genomename}.fna

	cp ${genomename}/${genomename}.gbk $OUTDIR/antismash
	cd $OUTDIR
	rm -f $genome
done

cd $OUTDIR/prokka_annotation
ls *.fna | rev | cut -f2- -d'.' | rev > $OUTDIR/genomes.txt

cd $OUTDIR/antismash

source activate antismash

for f in *.gbk
do
	antismash $f \
	--genefinding-tool none \
	--cb-subclusters --cb-knownclusters \
	--asf --cpus 6
	rm -f $f
done
conda deactivate


