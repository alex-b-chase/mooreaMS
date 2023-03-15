#!/bin/bash

#### get coverage stats for each MAG/bin for proxy for abundance
#### ideally want to correlate to MS data

#### can just do the dRep MAGs for convenience

BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly
MAGDIR=$BASEDIR/bowtie/combinedMAGs/testDrep/dereplicated_genomes

rm -rf $BASEDIR/tempcov
mkdir -p $BASEDIR/tempcov

cp $MAGDIR/*.fa $BASEDIR/tempcov
cp $BASEDIR/bowtie/sampleIDs.txt $BASEDIR/tempcov

cd $BASEDIR/tempcov

for f in *.fa
do
	binID=$(echo $f | sed 's/.bin./.B/g' | sed 's/.fa//g')
	while read mgID 
	do

		echo "#!/bin/bash
#PBS -q condo
#PBS -N ${binID}-${mgID}
#PBS -l nodes=1:ppn=4
#PBS -l walltime=07:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4

BASEDIR=$BASEDIR
OUTDIR=\$BASEDIR/tempcov

binID=${binID}
mgID=${mgID}

cd \$OUTDIR

bbmap.sh \\
in1=\$BASEDIR/\${mgID}.filter.clean.R1.fq.gz \\
in2=\$BASEDIR/\${mgID}.filter.clean.R2.fq.gz \\
ref=\$OUTDIR/$f \\
nodisk covstats=\${binID}.\${mgID}.covstats.txt

	" > $BASEDIR/tempcov/${binID}_${mgID}.sh
	# qsub $BASEDIR/tempcov/${binID}_${mgID}.sh

	done < sampleIDs.txt
done


#### after go through and collect data $ cat bbmap_cov.e29806043 | grep "Average coverage:" | cut -f2
