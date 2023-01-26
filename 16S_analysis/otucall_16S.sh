#!/bin/bash
#PBS -q home-sio
#PBS -N qiime2_MO18
#PBS -l nodes=1:ppn=36:mem768
#PBS -l walltime=18:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4

conda activate qiime2-2021.8

BASEDIR=/oasis/tscc/scratch/abchase/moorea/16S_analysis

### create 16S database and train the classifier
### only need to do once!
DBDIR=/oasis/tscc/scratch/abchase/refDB/16S_refDB

# cd $DBDIR
# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path 190206_ssu_r86.1_20180911.clean.fna \
#   --output-path GTDB_bac-arc_ssu_r86.qza

# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path ssu_r86.1_20180911.clean.qiime2.tax \
#   --output-path ref-taxonomy.qza

# ### trim database to the correct variable region of the 16S rRNA gene
# ### we used the 515F primer

# qiime feature-classifier extract-reads \
#   --i-sequences GTDB_bac-arc_ssu_r86.qza \
#   --p-f-primer GTGCCAGCMGCCGCGGTAA \
#   --p-r-primer GGACTACHVGGGTWTCTAAT \
#   --p-min-length 100 \
#   --p-max-length 400 \
#   --o-reads ref-seqs_515f806r.qza

# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads ref-seqs_515f806r.qza \
#   --i-reference-taxonomy ref-taxonomy.qza \
#   --o-classifier classifier_15f806r.qza

CLASSDB=$DBDIR/classifier_15f806r.qza


## DADA2 is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data
## “OTUs” resulting from DADA2 from unique sequences, these are the equivalent of 100% OTUs
cd $BASEDIR

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 9 \
  --p-trim-left-r 5 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table table-dada2.qza \
  --p-n-threads 36 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza

qiime metadata tabulate \
  --m-input-file denoising-stats-dada2.qza \
  --o-visualization stats-dada2.qzv

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file MO18-metadata.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

## make a phylogeny with repseqs
## need for some diversity metrics
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs-dada2.qza \
  --o-masked-alignment masked-aligned-rep-seqs-dada2.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

## assign taxonomy to the sequences
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSDB \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


conda deactivate
