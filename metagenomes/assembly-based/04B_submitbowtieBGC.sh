#!/bin/bash

BASEDIR=/dfs5/bio/abchase/moorea/idba
OUTDIR=$BASEDIR/bowtieBGCs

numthreads=8

cd $BASEDIR

while read MGID
do
	newsample=$(echo $MGID | sed 's/O18_//g')

	echo "#!/bin/bash
#$ -N ${newsample}_bowtieBGC
#$ -m a
#$ -q bio,mic
#$ -pe openmp ${numthreads}

module load bowtie2/2.2.7
module load samtools/1.3

BASEDIR=$BASEDIR
OUTDIR=\$BASEDIR/bowtieBGCs

BIN=/data/users/abchase/bin

MGID=$MGID
ASSEMBLY=BGCs_totalcontigs.fasta

cd \$OUTDIR

bowtie2 -x \$ASSEMBLY --very-fast-local -k 1 -t -p ${numthreads} --reorder --mm \\
-U <(\$BIN/shuffleSequences_fastx.pl 4 <(zcat \$BASEDIR/\${MGID}.fix.R1.fq.gz) <(zcat \$BASEDIR/\${MGID}.fix.R2.fq.gz)) \\
| samtools view -S -b -T \$ASSEMBLY - > \$OUTDIR/\$MGID.bowtie2.bam

	" > $OUTDIR/$MGID.bowtieBGC.sh
	# qsub $OUTDIR/$MGID.bowtieBGC.sh

done < $BASEDIR/sampleIDs.txt

# $BIN/gc_cov_annotate.pl \
# --blasttaxid rand.nt.1e-5.megablast \
# --assembly $ASSEMBLY --bam *.bam --out bowtie_blobplot.txt \
# --taxdump $TAXDUMP --taxlist genus family phylum


