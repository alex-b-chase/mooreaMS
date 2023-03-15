#!/bin/bash

BASEDIR=/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies

# rm -rf $BASEDIR/bigscape
# mkdir -p $BASEDIR/bigscape
# mkdir -p $BASEDIR/bigscape/input_files

# cd $BASEDIR/idba_output
# ls *.fna | cut -f1 -d'.' > $BASEDIR/genomes.txt

# cd $BASEDIR
# while read genome
# do
# 	echo "processing $genome ..."
# 	cd $BASEDIR/antismash/$genome 

# 	for f in *.region*.gbk 
# 	do
# 		cp $f $BASEDIR/bigscape/input_files/
# 	done

# done < genomes.txt


#################################################################
### run BiG-SCAPE and predict gene clusters families (GCFs)
#################################################################
cd $BASEDIR/bigscape

source activate bigscape

BIGDIR=/Users/alexchase/software/BiG-SCAPE
MIBIGDIR=/Volumes/JensenLabMGs/referenceDB/mibig_gbk_2.0

#### first get the really close hits to each other and MIBIG
python $BIGDIR/bigscape.py \
-i input_files \
-o bigscape_out \
--pfam_dir $BIGDIR \
-c 8 --mibig --cutoffs 0.8 \
--clans-off \
--mix --no_classify

#### need to copy over high hits from MIBIG
#### also add lobophorin for good measure
# cp $MIBIGDIR/BGC0001183.gbk $BASEDIR/bigscape/input_files/MB.region.lobophorinA.gbk
# cp $MIBIGDIR/BGC0001004.gbk $BASEDIR/bigscape/input_files/MB.region.lobophorinB.gbk



# python $BIGDIR/bigscape.py \
# -i input_files \
# -o bigscape_out \
# --pfam_dir $BIGDIR \
# -c 8 --cutoffs 0.8 \
# --clans-off \
# --mix --no_classify

conda deactivate 

