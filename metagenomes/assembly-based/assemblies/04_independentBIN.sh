#!/bin/bash

### Assemble each sample separately using your favorite assembler
### Bin each assembly (and co-assembly) separately using your favorite binner
### Pull the bins from all assemblies together and run dRep on them

BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly
OUTDIR=$BASEDIR/bowtie

numthreads=8

rm -rf $OUTDIR/combinedMAGs
mkdir -p $OUTDIR/combinedMAGs

while read MGID 
do

	echo "#!/bin/bash
#PBS -q home-sio
#PBS -N ${MGID}_bowtie2
#PBS -l nodes=1:ppn=8:mem128
#PBS -l walltime=47:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4
module load bowtie2/2.3.4.3
module load samtools/1.9
module load hmmer/3.2.1 

BASEDIR=${BASEDIR}
OUTDIR=${OUTDIR}

MGID=${MGID}
ASSEMBLY=\$BASEDIR/contigs/\${MGID}.prod.fa

cd \$OUTDIR

bowtie2-build \$ASSEMBLY \$ASSEMBLY

while read sampleID 
do

	bowtie2 -x \$ASSEMBLY --sensitive -p ${numthreads} --maxins 800 \\
	-1 \$BASEDIR/\${sampleID}.filter.clean.R1.fq.gz -2 \$BASEDIR/\${sampleID}.filter.clean.R2.fq.gz \\
	| samtools view -S -b -T \$ASSEMBLY - > \$OUTDIR/\${MGID}.\${sampleID}.bowtie2.bam

	samtools sort \$OUTDIR/\${MGID}.\${sampleID}.bowtie2.bam \\
	--threads ${numthreads} > \$OUTDIR/\${MGID}.\${sampleID}.sorted.bam

	rm -f \$OUTDIR/\${MGID}.\${sampleID}.bowtie2.bam

done < \$OUTDIR/sampleIDs.txt

module load gnu
export CC=/opt/gnu/gcc/bin/gcc
export CXX=/opt/gnu/gcc/bin/g++

MINLEN=5000

jgi_summarize_bam_contig_depths \\
--outputDepth \$OUTDIR/\${MGID}.metabat_depth.txt \\
\$OUTDIR/\${MGID}.*.sorted.bam

rm -rf \${MGID}.metabat/

metabat2 -i \$ASSEMBLY \\
-a \$OUTDIR/\$MGID.metabat_depth.txt -o \${MGID}.metabat/\${MGID}.bin \\
-m \$MINLEN -t ${numthreads} -v

cp \${MGID}.metabat/*.fa \$OUTDIR/combinedMAGs

rm -f \$OUTDIR/\${MGID}.*.sorted.bam

	" > $OUTDIR/$MGID.bowtie.sh
	qsub $OUTDIR/$MGID.bowtie.sh
	sleep 10s

done < $OUTDIR/sampleIDs.txt




