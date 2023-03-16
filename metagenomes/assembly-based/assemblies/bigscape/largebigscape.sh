#!/bin/bash

REFDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/bigscape 

# a lot of the BGCs that were identified were small and fragmented
# let's take larger BGCs (>10 genes) and recreate the BiG-SCAPE plot
### FOR THE PUBLICATION, USE ALL BGCS FOR NETWORK

# mkdir -p $REFDIR/input_files10
# cd $REFDIR/input_files

# mincount=10

# for f in *.gbk
# do
# 	BGCID=${f%.gbk}

# 	genbank_to_fasta.py -i $f -s aa -o ${BGCID}.faa 

# 	genecount=$(grep ">" ${BGCID}.faa | wc -l)

# 	if [ "$genecount" -gt "$mincount" ]
# 	then
# 		cp $f $REFDIR/input_files10
# 	else
# 		:
# 	fi

# 	rm -f ${BGCID}.faa 
# done

cd $REFDIR

source activate bigscape

BIGDIR=/opt/miniconda3/envs/BiG-SCAPE
PFAMDIR=/Users/alexchase/software/PfamScan

MIBIGDIR=/Volumes/JensenLabMGs/referenceDB/mibig_gbk_2.0

# #### first get the really close hits to each other and MIBIG
# #### then run with mibig hits from 0.5 threshold
# #### also add lobophorin for good measure
rm -f $REFDIR/input_files/MB.*
# cp $MIBIGDIR/BGC0001183.gbk $REFDIR/input_files10/MB.region.lobophorinA.gbk
# cp $MIBIGDIR/BGC0001004.gbk $REFDIR/input_files10/MB.region.lobophorinB.gbk

# cd $REFDIR/bigscape_out/network_files/2021-06-18_14-03-55_hybrids_glocal

# cut -f1 mix/mix_clustering_c0.50.tsv | grep "^BGC" | cut -f1 -d'.' > $REFDIR/mibighits.txt

mappingfile=$REFDIR/bigscape_out/network_files/2021-06-18_16-02-21_hybrids_glocal/Network_Annotations_Full.tsv

while read line
do
	BGCname=$(grep "$line" $mappingfile | cut -f3 | rev | cut -f4- -d' ' | rev | sed 's/ //g' | tr -d "[ -%,;\(\):=\.\\\[]\"\']--//")
	cp $MIBIGDIR/${line}.gbk $REFDIR/input_files/MB.region.${line}.gbk
done < $REFDIR/mibighits.txt

cd $REFDIR

python $BIGDIR/bigscape.py \
-i input_files \
-o bigscape_out \
--pfam_dir $PFAMDIR \
-c 8 --cutoffs 0.6 \
--clans-off \
--mix --no_classify

conda deactivate 
