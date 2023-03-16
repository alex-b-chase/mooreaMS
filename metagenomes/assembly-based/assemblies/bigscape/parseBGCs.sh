#!/bin/bash

BASE=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies
BASEDIR=$BASE/bigscape
BIGSCAPEDIR=$BASEDIR/input_files

# info to make the attribute excel sheet for network in GEPHI or Cytoscape
cd $BASEDIR
rm -f $BASEDIR/contigedges.txt
rm -f $BASEDIR/temp.txt
rm -f $BASEDIR/mibighits.txt
## first get the mibig reference clusters
cut -f1 -d',' bgc_cluster_network.csv | tr -d '"' | grep -v "Source" > temp.txt
cut -f2 -d',' bgc_cluster_network.csv | tr -d '"' | grep -v "Target" >> temp.txt

cat temp.txt | sort -u > network.nodes.uniq.txt 
cat network.nodes.uniq.txt | grep "^MB"  | cut -f3 -d'.' > mibig.temp.txt
rm -f temp.txt

while read BGC
do 
	grep $BGC $BASEDIR/bigscape_out/network_files/2022-10-20_12-08-51_hybrids_glocal/Network_Annotations_Full.tsv >> mibighits.txt
done < mibig.temp.txt 
rm -f mibig.temp.txt

cut -f1,3 mibighits.txt | sed 's/ biosynthetic gene cluster//g' | sed 's/ //g' | awk '{sub(/\.1.*$/,"",$1)}1' OFS='\t' > mibig_BGC_innetwork.reduce.txt 

# get more info on whether the BGCs are in complete or on contig edges
# also check whether the BGCs are in med-high quality MAGs
cd $BASE/binned_MAGs/binned_genomes
cat *.fa > total.contigs.fasta 
print-fasta-id.py total.contigs.fasta 
mv total.contigs_ids.txt $BASEDIR/ 
rm -f total.contigs.fasta 

mags2contigs=$BASEDIR/total.contigs_ids.txt

cd $BIGSCAPEDIR
echo -e "Id\tBGCtype\tComplete\tMAGhit\tMAG" > $BASEDIR/BGCinformation.txt
for f in *.gbk
do 
	contigedge=$(grep 'contig_edge="False"' $f | wc -l)
	BGCtype=$(grep "/product=" $f | head -n1 | sed 's/^ *//g' | sed 's;/product=;;g' | tr -d '"')
	BGCname=${f%.gbk}

	contigname=$(echo $BGCname | cut -f1-2 -d'.')
	MAGcontained=$(grep -w "${contigname}" $mags2contigs )

	if [ $contigedge -gt 1 ]
	then
		contighit=$(echo "YES")
	else 
		contighit=$(echo "NO")
	fi

	if [ -n "$MAGcontained" ]
	then
		MAGhit=$(echo "YES")
		for binMAG in $BASE/binned_MAGs/binned_genomes/*.fa 
		do
			result=$(grep -w ">${contigname}" $binMAG )
			if [ -n "$result" ]; then
				whichMAG=$(echo $binMAG | rev | cut -f1 -d'/' | rev | cut -f1-3 -d'.')
				break
			else
				whichMAG=$(echo "NA")
			fi
		done
	else 
		MAGhit=$(echo "NO")
		whichMAG=$(echo "NA")
	fi
	echo -e "${BGCname}\t${BGCtype}\t${contighit}\t${MAGhit}\t${whichMAG}" >> $BASEDIR/BGCinformation.txt
done

rm -f total.contigs_ids.txt


