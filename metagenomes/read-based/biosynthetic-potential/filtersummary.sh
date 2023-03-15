#!/bin/bash

### make reference databases from NaPDoS v2 files
### downloaded from http://npdomainseeker.sdsc.edu/napdos2/

REFDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/gene-centric-approach/biosynthetic_potential
KSDB=$REFDIR/napdosV2_KS
ConDB=$REFDIR/napdosV2_C

# for f in *.faa
# do
# 	clustalo -i $f -o ${f%.faa}.aln

# 	# raxml -s ${f%.faa}.aln -m PROTGAMMAWAGF -n ${f%.faa} -x 100 -# 100 -p 4321 -f a -T 6

# 	hmmbuild ${f%.faa}.hmm ${f%.faa}.aln

# done

### now HMMer for the KS and C domains

cd $REFDIR/blat_output

rm -f *.hmm.faa
rm -f *.temp.*

rm -rf summaryoutput/
mkdir -p summaryoutput
# rm -rf tempdir
mkdir -p tempdir/
mkdir -p summaryoutput

ls *.faa | cut -f1 -d'.' | cut -f1-2 -d'_' | sort -u | grep "^M" > sampleIDs.txt

echo -e "sampleID\tBGC\ttotalBGChits.blat\ttotalGChits.hmmer" > $REFDIR/summaryBP.txt


while read sampleID
do
	metagenome=$(echo $sampleID | sed 's/MO18/MO/g')

	# first filter out C domains
	for BGC in C KS
	do

		cat ${sampleID}*.blat.${BGC}markers.faa | \
		sed -e "s/^>/>${metagenome}_/" | sed '/^>/ s/ .*//' > ${metagenome}.${BGC}.faa

		mv ${sampleID}*.blat.${BGC}markers.faa tempdir

		blathits=$(grep ">" ${metagenome}.${BGC}.faa | wc -l)

		hmmsearch --tblout ${metagenome}.${BGC}.hmm.txt -E 1e-10 \
		--cpu 4 $REFDIR/napdosV2_${BGC}domains.hmm ${metagenome}.${BGC}.faa > ${metagenome}.${BGC}.log.txt

		cat ${metagenome}.${BGC}.hmm.txt | cut -f1 -d' ' | sort | uniq > ${metagenome}.${BGC}.temp.txt

		# subset new HMMer filtered reads
		filterbyname.sh \
		in=${metagenome}.${BGC}.faa \
		out=${metagenome}.${BGC}.hmm.faa \
		names=${metagenome}.${BGC}.temp.txt \
		ow=t include=t tossjunk 2>/dev/null

		cat ${metagenome}.${BGC}.hmm.faa | tr -d '*' > ${metagenome}.${BGC}.hmm.temp.faa
		hmmerhits=$(grep ">" ${metagenome}.${BGC}.hmm.temp.faa | wc -l)

		echo -e "${metagenome}\t${BGC}\t${blathits}\t${hmmerhits}" >> $REFDIR/summaryBP.txt

		mv ${metagenome}.${BGC}.faa summaryoutput
		rm -f ${metagenome}.${BGC}.temp.txt
		rm -f ${metagenome}.${BGC}.hmm.txt
		rm -f ${metagenome}.${BGC}.log.txt
		rm -f ${metagenome}.${BGC}.hmm.faa

	done

done < sampleIDs.txt

cat *.C.hmm.temp.faa > total.Cdomains.faa
cat *.KS.hmm.temp.faa > total.KSdomains.faa

rm -f *.C.hmm.temp.faa
rm -f *.KS.hmm.temp.faa

#### run the total.*domains.faa through NaPDoSv2

source activate biotools

diamond blastp \
--query total.Cdomains.faa \
--db $ConDB --min-orf 50 --evalue 0.00001 | \
sort -k1,1 -k2,2 -k12,12nr | sort -u -k1,1 > total.Cdomains.topblast.txt

diamond blastp \
--query total.KSdomains.faa \
--db $KSDB --min-orf 50 --evalue 0.00001 | \
sort -k1,1 -k2,2 -k12,12nr | sort -u -k1,1 > total.KSdomains.topblast.txt

conda deactivate

cut -f1-2 total.Cdomains.topblast.txt > $REFDIR/napdos_output/cdomains_TOTAL.txt
cut -f1-2 total.KSdomains.topblast.txt > $REFDIR/napdos_output/ksdomains_TOTAL.txt

