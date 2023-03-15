#!/bin/bash

BASEDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing
KSDB=$BASEDIR/gene-centric-approach/biosynthetic_potential/napdosV2_KS
ConDB=$BASEDIR/gene-centric-approach/biosynthetic_potential/napdosV2_C

MGDIR=$BASEDIR/MAGs/assemblies

cd $MGDIR

echo -e "sampleID\toldAssembliesKS\tKSdomainsL2kbp\tKSdomainsL5kbp" > $MGDIR/napdoshits.txt

source activate biotools

## only need to do once
### diamond makedb --in napdosV2_KSdomains.faa --db napdosV2_KS 

while read sampleID
do

	diamond blastx \
	--query $MGDIR/${sampleID}.contigs.L2kbp.fna \
	--db $KSDB --min-orf 200 --evalue 0.00001 | \
	sort -k1,1 -k2,2 -k12,12nr | sort -u -k1,1 > ${sampleID}.topblast.L2kbp.txt

	L2kbp=$(cat ${sampleID}.topblast.L2kbp.txt | wc -l)

	diamond blastx \
	--query $MGDIR/antismash/${sampleID}.prod.fa \
	--db $KSDB --min-orf 200 --evalue 0.00001 | \
	sort -k1,1 -k2,2 -k12,12nr | sort -u -k1,1 > ${sampleID}.topblast.L5kbp.txt

	L5kbp=$(cat ${sampleID}.topblast.L5kbp.txt | wc -l)

	diamond blastx \
	--query $BASEDIR/MAGs/assemblies-OLD/fixed_assemblies/ncbi_fix/${sampleID}.contigs.L1kbp.fna \
	--db $KSDB --min-orf 200 --evalue 0.00001 | \
	sort -k1,1 -k2,2 -k12,12nr | sort -u -k1,1 > ${sampleID}.topblast.old.txt

	oldass=$(cat ${sampleID}.topblast.old.txt | wc -l)

	echo -e "${sampleID}\t${oldass}\t${L2kbp}\t${L5kbp}" >> $MGDIR/napdoshits.txt

done < genomes.txt

cat *.topblast.L2kbp.txt > napdos.KSdomains.L2kbp.txt
cat *.topblast.L5kbp.txt > napdos.KSdomains.L5kbp.txt

rm -f *.topblast.L2kbp.txt
rm -f *.topblast.L5kbp.txt

