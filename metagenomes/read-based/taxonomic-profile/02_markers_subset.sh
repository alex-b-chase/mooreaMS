#!/bin/bash

# take the output from BLAT results and subset into marker gene specific files
# for some reason, (think the memory alottment), the HPC runs this SUPER slowly
# use locally and then do the HMMER and PPLACER analyses on the HPC

BASEDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/gene-centric-approach/taxonomic-profile
OUTDIR=$BASEDIR/blat_output

# subset the BLAT results from the MG libraries

cd $OUTDIR
cp /Volumes/JensenLabMGs/referenceDB/comm_markers/marker_genes.txt $BASEDIR/

rm -f *.blat.faa
rm -f *.temp.*

rm -rf summaryoutput
mkdir -p summaryoutput
# rm -rf tempdir
mkdir -p tempdir
mkdir -p summaryoutput

ls *.txt | cut -f1 -d'.' | cut -f1-2 -d'_' | sort -u | grep "^M" > sampleIDs.txt

while read sampleID
do

	metagenome=$(echo $sampleID | sed 's/MO18/MO/g')

	cat ${sampleID}*.markers.faa | \
	sed -e "s/^>/>${metagenome}_/" | sed '/^>/ s/ .*//' > ${metagenome}.totalmarkers.temp.faa
	cat ${sampleID}*.blat.total.txt | \
	sed -e "s/^/${metagenome}_/" | cut -f1-2 > ${metagenome}.blat.temp.txt

	cat ${sampleID}*.markers.faa > summaryoutput/${metagenome}.markers.faa
	cat ${sampleID}*.blat.total.txt > summaryoutput/${metagenome}.blat.total.txt
	mv ${sampleID}*.markers.faa tempdir
	mv ${sampleID}*.blat.total.txt tempdir

done < sampleIDs.txt

cat *.temp.faa > total.markers.faa 
cat *.blat.temp.txt > total.blatmarkers.txt 

rm -f *.temp.faa
rm -f *.blat.temp.txt

# now subset by each marker gene

while read protein
do

	echo "Processing ${protein}..."

	
	if [ ! -f "${protein}.blat.faa" ]
	then

		grep -w "${protein}p" total.blatmarkers.txt | cut -f1 | sort -u | sed 's/[ ]*$//' > ${protein}.temp.txt

		# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
		filterbyname.sh \
		in=total.markers.faa  \
		out=${protein}.temp.faa \
		names=${protein}.temp.txt ow=t include=t tossjunk 2>/dev/null
		rm ${protein}.temp.txt

		# rename each read to index each sequence
		cat ${protein}.temp.faa | awk '/^>/ {$0=NR"_"$0} 1' | sed 's/ //g; s/>//g' | awk '/^[1-9]/ {$0=">"$0} 1' > ${protein}.blat.faa
		rm ${protein}.temp.faa

	else 
		continue 
	fi

done < $BASEDIR/marker_genes.txt





