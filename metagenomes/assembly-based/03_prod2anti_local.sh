#!/bin/bash

BASEDIR=/Volumes/JensenLabMGs/alex_alyssa/moorea2020/metagenome_processing/MAGs/assemblies

# ran IDBA assembly - now need to run through antiSMASH v5 for BGC prediction
# probably want to filter to <5kbp contigs

mkdir -p $BASEDIR/antismash

cd $BASEDIR

# first rename the contig files so antiSMASH likes them
export PATH="/anaconda2/bin:$PATH"
eval "$(conda shell.bash hook)"

for f in *.L2kbp.fna
do
	sampleID=${f%.contigs.L2kbp.fna}
	newsampleID=$(echo $sampleID | sed 's/MO18_/M/g')
	echo "${sampleID}"

	if [ -d "$BASEDIR/antismash/${sampleID}" ]; then
		echo "done with ${sampleID}"
		continue
	else
		# rm -f $BASEDIR/antismash/${sampleID}.prod.gff
		# rm -f $BASEDIR/antismash/${sampleID}.prod.fa

		n50_calc.py $f 

		rename.sh in=$f out=${sampleID}.temp.fa prefix=${newsampleID}.c

		bbduk.sh in=${sampleID}.temp.fa \
		out=$BASEDIR/antismash/${sampleID}.prod.fa \
		minlen=5000 ow=t

		rm -f ${sampleID}.temp.fa
		n50_calc.py $BASEDIR/antismash/${sampleID}.prod.fa

		prodigal -i $BASEDIR/antismash/${sampleID}.prod.fa \
		-o $BASEDIR/antismash/${sampleID}.prod.gff \
		-p meta -f gff 

		conda activate antismash
		antismash \
		--genefinding-gff3 $BASEDIR/antismash/${sampleID}.prod.gff \
		--output-dir $BASEDIR/antismash/${sampleID} \
		$BASEDIR/antismash/${sampleID}.prod.fa
		conda deactivate
	fi

done

# rm -rf $BASEDIR/bigscape
# mkdir -p $BASEDIR/bigscape
# mkdir -p $BASEDIR/bigscape/input_files

# ls *.fna | cut -f1 -d'.' > genomes.txt

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

conda activate bigscape

BIGDIR=/anaconda2/envs/bigscape/BiG-SCAPE
PFAMDIR=/Users/alexchase/software/PfamScan

python $BIGDIR/bigscape.py \
-i input_files \
-o bigscape_out \
--pfam_dir $PFAMDIR \
-c 4 --mibig --cutoffs 0.8 \
--clan_cutoff 0.8 0.9 \
--mix --no_classify

conda deactivate 

