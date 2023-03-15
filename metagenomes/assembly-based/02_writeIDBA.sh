#!/bin/bash

THREAD=16
BASE=/oasis/tscc/scratch/abchase/moorea
BASEDIR=$BASE/assembly

cd $BASEDIR

rm -rf $BASEDIR/idba
mkdir -p $BASEDIR/idba

while read sampleID
do

	# need to fix the reads since the first pass did not remove all adapter sequences...
	# SUPER FRUSTRATINGGGSDG

	# for replicate in A B C
	# do

	# 	bbduk.sh qtrim=rl trimq=10 threads=2 \
	# 	minlen=25 ktrim=r k=25 ref=/data/users/abchase/adapters.bbmap.fa hdist=1 \
	# 	in1=${sampleID}_${replicate}.filter.clean.R1.fq.gz in2=${sampleID}_${replicate}.filter.clean.R2.fq.gz \
	# 	out1=$sampleID.fix.fq.gz out2=$sampleID.fix2.fq.gz tbo

	# 	repair.sh in=$sampleID.fix.fq.gz in2=$sampleID.fix2.fq.gz \
	# 	out=${sampleID}_${replicate}.fix.R1.fq.gz out2=${sampleID}_${replicate}.fix.R2.fq.gz

	# 	rm -f $sampleID.fix.fq.gz
	# 	rm -f $sampleID.fix2.fq.gz
	# done

	# submit the HPC jobs for each sampleID
	echo "#!/bin/bash
#PBS -q home-sio
#PBS -N ${sampleID}_IDBA
#PBS -l nodes=1:ppn=${THREAD}:mem128
#PBS -l walltime=23:59:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

REF=${sampleID}
REFBASE=${BASEDIR}
OUTDIR=\$REFBASE/idba

## we have technical replicates of the sampleIDs in 3 different copies
## combine these 3 technical replicates and assemble from there from initial samples 
##cat \${REF}_A.fix.R1.fq.gz \${REF}_B.fix.R1.fq.gz \${REF}_C.fix.R1.fq.gz > \$OUTDIR/\${REF}.fix.R1.fq.gz 
##cat \${REF}_A.fix.R2.fq.gz \${REF}_B.fix.R2.fq.gz \${REF}_C.fix.R2.fq.gz > \$OUTDIR/\${REF}.fix.R2.fq.gz 

## having BIG memory issues with running IDBA assembler
## not too surprising but need to digitally normalize the libraries
## i.e., many assemblers perform poorly in the presence of too much data, 
## and data with irregular coverage, such as MDA-amplified single cells or metagenomes

cd \$REFBASE
repair.sh in=\${REF}.filter.clean.R1.fq.gz in2=\${REF}.filter.clean.R2.fq.gz \\
out=\$OUTDIR/\${REF}.fix.R1.fq.gz out2=\$OUTDIR/\${REF}.fix.R2.fq.gz ow=t

cd \$OUTDIR

bbnorm.sh \\
in=\${REF}.fix.R1.fq.gz in2=\${REF}.fix.R2.fq.gz \\
target=40 mindepth=5 -Xmx32G ow=t \\
out=\${REF}.norm.R1.fq.gz out2=\${REF}.norm.R2.fq.gz

rm -f \${REF}.fix.R1.fq.gz
rm -f \${REF}.fix.R2.fq.gz

cd \$OUTDIR

conda activate idba
## IDBA needs the reads to be in fastA format 

fq2fa --merge <(zcat \${REF}.norm.R1.fq.gz) <(zcat \${REF}.norm.R2.fq.gz) \${REF}.fas

rm -rf \$OUTDIR/\${REF}/
rm -f \${REF}.norm.R1.fq.gz
rm -f \${REF}.norm.R2.fq.gz

idba_ud -r \${REF}.fas --pre_correction \\
--mink 30 --maxk 200 --step 10 --num_threads ${THREAD} \\
--min_contig 2000 --out \${REF}

conda deactivate

## extra precaution to gets reads into suitable format for binning steps
bbduk.sh in=\${REF}/contig.fa out=\$OUTDIR/\${REF}.contigs.L2kbp.temp.fna minlen=2000 ow=t

rename.sh in=\$OUTDIR/\${REF}.contigs.L2kbp.temp.fna \\
out=\$OUTDIR/\${REF}.contigs.L2kbp.fna prefix=\${REF} addprefix=t ow=t

rm -f \$OUTDIR/\${REF}.contigs.L2kbp.temp.fna

##rm -rf \$OUTDIR/\${REF}/

	" > $BASEDIR/idba/$sampleID.idba.sh

	# qsub $BASEDIR/idba/$sampleID.idba.sh

done < $BASE/sampleIDs.txt

