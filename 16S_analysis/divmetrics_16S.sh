#!/bin/bash
#PBS -q home-sio
#PBS -N qiime2_divmet
#PBS -l nodes=1:ppn=8:mem768
#PBS -l walltime=02:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4

conda activate qiime2-2021.8

BASEDIR=/oasis/tscc/scratch/abchase/moorea/16S_analysis
DBDIR=/oasis/tscc/scratch/abchase/refDB/16S_refDB
CLASSDB=$DBDIR/classifier_15f806r.qza

cd $BASEDIR

rm -f taxa-bar-plots.qzv
rm -f filtered-table-temp.qza
rm -f filtered-table-dada2.qza
rm -f filtered-table-dada2.qzv
rm -rf core-metrics-results-all/

## need to get rid of chloroplasts! might as well do mitochondria & eukaryota too although not usually an issue
## also remove singleton ESVs
qiime taxa filter-table \
--i-table table-dada2.qza \
--i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast,eukaryota \
--o-filtered-table filtered-table-temp.qza

#filter singletons - maybe should not do? according to Jen, should keep before rarefaction
# easy to remove later, remove for now
# qiime feature-table filter-features \
# --i-table filtered-table-temp.qza \
# --p-min-frequency 2 \
# --o-filtered-table filtered-table-allsamps.qza

# ### remove transect samples from data to see if that helps ordination - doesn't help!
# qiime feature-table filter-samples \
# --i-table filtered-table-allsamps.qza \
# --m-metadata-file MO18-metaQuad.txt \
# --o-filtered-table filtered-table-dada2.qza

### remove site 6 samples from data shouldn't have been sampled anyways....
qiime feature-table filter-samples \
--i-table filtered-table-allsamps.qza \
--m-metadata-file MO18-metaQuad.txt \
--o-filtered-table filtered-table-dada2.qza

rm -f filtered-table-temp.qza

#visualize
qiime feature-table summarize \
--i-table filtered-table-allsamps.qza \
--o-visualization filtered-table-allsamps.qzv \
--m-sample-metadata-file MO18-metadata.txt

qiime tools export \
--input-path filtered-table-allsamps.qza \
--output-path OTUtable-all

biom convert -i OTUtable-all/feature-table.biom -o feature-table-all.tsv --to-tsv

## simple barplot across samples
qiime taxa barplot \
  --i-table filtered-table-allsamps.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file MO18-metadata.txt \
  --o-visualization taxa-bar-plots-all.qzv

## alpha diversity 
## Shannon’s diversity index (a quantitative measure of community richness)

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table-allsamps.qza \
  --p-sampling-depth 17559 \
  --m-metadata-file MO18-metadata.txt \
  --output-dir core-metrics-results-all

## explore the microbial composition of the samples in the context of the sample metadata
## associations between categorical metadata columns and alpha diversity data

## Faith Phylogenetic Diversity (a measure of community richness)
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-all/faith_pd_vector.qza \
  --m-metadata-file MO18-metadata.txt \
  --o-visualization core-metrics-results-all/faith-pd-group-significance.qzv

## evenness metrics
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-all/evenness_vector.qza \
  --m-metadata-file MO18-metadata.txt \
  --o-visualization core-metrics-results-all/evenness-group-significance.qzv

## beta diversity
## sample composition in the context of categorical metadata using PERMANOVA

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-all/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MO18-metadata.txt \
  --m-metadata-column site_number \
  --o-visualization core-metrics-results-all/weighted-unifrac-site-significance.qzv \
  --p-pairwise

## plot ordination 
## use the Emperor tool to explore principal coordinates (PCoA) plots in the context of sample metadata

qiime emperor plot \
  --i-pcoa core-metrics-results-all/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file MO18-metadata.txt \
  --o-visualization core-metrics-results-all/unweighted-unifrac-emperorplot.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results-all/bray_curtis_pcoa_results.qza \
  --m-metadata-file MO18-metadata.txt \
  --o-visualization core-metrics-results-all/bray-curtis-emperorplot.qzv

conda deactivate
