#!/bin/bash

THREAD=2
BASEDIR=/dfs5/bio/abchase/moorea
OUTDIR=$BASEDIR/idba

cd $OUTDIR

while read sampleID
do
	jobname=$(echo $sampleID | sed 's/O18_//g')

	echo "#!/bin/bash
#$ -N ${jobname}.combine
#$ -m a
#$ -j y
#$ -q bio
#$ -pe openmp ${THREAD}

module load BBMap/37.50

REF=${sampleID}
REFBASE=${BASEDIR}
OUTDIR=\$REFBASE/idba

cd \$REFBASE

## we have technical replicates of the sampleIDs in 3 different copies
## combine these 3 technical replicates and assemble from there

## having BIG memory issues with running IDBA assembler
## not too surprising but need to digitally normalize the libraries
## i.e., many assemblers perform poorly in the presence of too much data, 
## and data with irregular coverage, such as MDA-amplified single cells or metagenomes

cat \${REF}_A.fix.R1.fq.gz \${REF}_B.fix.R1.fq.gz \${REF}_C.fix.R1.fq.gz > \$OUTDIR/\${REF}.fix.R1.fq.gz 
cat \${REF}_A.fix.R2.fq.gz \${REF}_B.fix.R2.fq.gz \${REF}_C.fix.R2.fq.gz > \$OUTDIR/\${REF}.fix.R2.fq.gz 

	" > $OUTDIR/$jobname.bowtie.sh


done < $BASEDIR/sampleIDs.txt



