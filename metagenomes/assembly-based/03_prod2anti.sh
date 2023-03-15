#!/bin/bash

BASE=/oasis/tscc/scratch/abchase/moorea/assembly
OUTDIR=$BASE/antismash

rm -rf $OUTDIR
mkdir -p $OUTDIR

cd $BASE/idba

for metagenome in *.contigs.L2kbp.fna
do

	genomename=${metagenome%.contigs.L2kbp.fna}

	echo "#!/bin/bash
#PBS -q condo
#PBS -N ${genomename}_prod2anti
#PBS -l nodes=1:ppn=8:mem128
#PBS -l walltime=7:59:00
#PBS -A jensen-group

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

BASEDIR=${BASE}
OUTDIR=${OUTDIR}
INPUTFILE=${metagenome}
sampleID=${genomename}

cd \${BASEDIR}/idba

rename.sh in=\$INPUTFILE \\
out=\${sampleID}.temp.fa prefix=\${sampleID}.c

bbduk.sh in=\${sampleID}.temp.fa \\
out=\$OUTDIR/\${sampleID}.prod.fa \\
minlen=5000 ow=t

rm -f \${sampleID}.temp.fa
cd \$OUTDIR

prodigal \\
-i \${sampleID}.prod.fa \\
-o \${sampleID}.prod.gff \\
-p meta -f gff 

conda activate antismash

export JAVA_OPTS=\"-Xmx2048m -XX:CompressedClassSpaceSize=256m\"
export _JAVA_OPTIONS=\"-Xmx2048m -XX:CompressedClassSpaceSize=256m\"

antismash \
--genefinding-gff3 \${sampleID}.prod.gff \\
--output-dir \${sampleID} \\
--cb-subclusters --cb-knownclusters \\
\${sampleID}.prod.fa

rm -f \${sampleID}.prod.gff

conda deactivate

tar -zcvf \${sampleID}.tar.gz ./\${sampleID}/

	" > $OUTDIR/$genomename.anti.sh

	qsub $OUTDIR/$genomename.anti.sh
	sleep 2s

done


