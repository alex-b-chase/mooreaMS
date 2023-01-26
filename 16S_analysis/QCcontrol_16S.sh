#!/bin/bash
#PBS -q home-sio
#PBS -N qiime2_MO18
#PBS -l nodes=1:ppn=4
#PBS -l walltime=03:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4

conda activate qiime2-2021.8

BASEDIR=/oasis/tscc/scratch/abchase/moorea/16S_analysis
OUTDIR=$BASEDIR/combinedseq

rm -rf $OUTDIR
mkdir -p $OUTDIR

## we have two sequencing runs for 16S data, might as well combine data
REFDIR1=$BASEDIR/samples
REFDIR2=$BASEDIR/samples2

cd $BASEDIR

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $BASEDIR/manifest.tsv

while read sampleID
do

	cat $REFDIR1/${sampleID}*_R1_*.fastq.gz \
	$REFDIR2/${sampleID}*_R1_*.fastq.gz > $OUTDIR/${sampleID}_R1_001.fastq.gz
	cat $REFDIR1/${sampleID}*_R2_*.fastq.gz \
	$REFDIR2/${sampleID}*_R2_*.fastq.gz > $OUTDIR/${sampleID}_R2_001.fastq.gz

	repair.sh \
	in=$OUTDIR/${sampleID}_R1_001.fastq.gz \
	in2=$OUTDIR/${sampleID}_R2_001.fastq.gz \
	out=$OUTDIR/${sampleID}_R1.fix.fastq.gz \
	out2=$OUTDIR/${sampleID}_R2.fix.fastq.gz \
	outs=$OUTDIR/${sampleID}.singleton.fastq.gz

	rm -f $OUTDIR/${sampleID}_R1_001.fastq.gz
	rm -f $OUTDIR/${sampleID}_R2_001.fastq.gz

	echo -e "${sampleID}\t${OUTDIR}/${sampleID}_R1.fix.fastq.gz\t${OUTDIR}/${sampleID}_R2.fix.fastq.gz" >> $BASEDIR/manifest.tsv

done < sampleIDs.txt

cd $BASEDIR

rm -f demux*

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

conda deactivate
