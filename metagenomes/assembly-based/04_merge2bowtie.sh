#!/bin/bash
#PBS -N cd2blast
#PBS -q home-sio
#PBS -l nodes=1:ppn=36:mem768
#PBS -l walltime=55:59:00

BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly
OUTDIR=$BASEDIR/bowtie

mkdir -p $OUTDIR

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4
module load bowtie2/2.3.4.3
module load samtools/1.9

numthreads=36

BIN=/home/abchase/bin
BLASTDB=/oasis/tscc/scratch/abchase/refDB/blastDB
TAXDUMP=/oasis/tscc/scratch/abchase/refDB/taxdump_ncbi

cd $BASEDIR/antismash

projectname=merged.assemblies
merged=${projectname}.fasta
merged_clustered=${projectname}.s99.fasta
afg_format=${projectname}.s99.afg

rm -f $OUTDIR/${projectname}.log

ls *.prod.fa | cut -f1 -d'.' > $OUTDIR/sampleIDs.txt

# combine the filtered assemblies together for dereplication
cat *.prod.fa > $OUTDIR/$merged

cd $OUTDIR
# # take the individual assemblies; identify and remove contigs contained in other assemblies
## needed to install newer version for the increased RAM power
CDHITEXEC=/home/abchase/software/cdhit-master

## need to load these for some reason....
module load gnu
export CC=/opt/gnu/gcc/bin/gcc
export CXX=/opt/gnu/gcc/bin/g++

$CDHITEXEC/cd-hit-est -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9 

### now ready to map to the assembled contigs
ASSEMBLY=$merged_clustered

## sample some random contigs for taxonomic hits - take too long to do all of them
$BIN/fastaqual_select.pl -f $ASSEMBLY -s r -n 10000 > random.${ASSEMBLY}

conda activate blast 

blastn -task megablast \
-query random.${ASSEMBLY} -db $BLASTDB/nt -evalue 1e-5 -max_target_seqs 1 \
-num_threads $numthreads -outfmt '6 qseqid staxids' -out $OUTDIR/rand.nt.1e-5.megablast

# index the assembly for mapping for the coverages
bowtie2-build $ASSEMBLY $ASSEMBLY

while read MGID 
do

	echo "#!/bin/bash
#PBS -q condo
#PBS -N ${MGID}_bowtie2
#PBS -l nodes=1:ppn=8:mem128
#PBS -l walltime=7:59:00
#PBS -A jensen-group

source ~/.bashrc
 
module load intel/2018.1.163
module load openmpi_ib/3.1.4
module load bowtie2/2.3.4.3
module load samtools/1.9

BASEDIR=${BASEDIR}
OUTDIR=\$BASEDIR/bowtie

BIN=${BIN}

MGID=${MGID}
ASSEMBLY=${ASSEMBLY}

cd \$OUTDIR

bowtie2 -x \$ASSEMBLY --very-fast-local -k 1 -t -p 8 --reorder --mm \\
-1 \$BASEDIR/\${MGID}.filter.clean.R1.fq.gz -2 \$BASEDIR/\${MGID}.filter.clean.R2.fq.gz \\
| samtools view -S -b -T \$ASSEMBLY - > \$OUTDIR/\$MGID.bowtie2.bam

	" > $OUTDIR/$MGID.bowtie.sh
	# qsub $OUTDIR/$MGID.bowtie.sh

done < $OUTDIR/sampleIDs.txt

conda deactivate


#### run this after all jobs are done
# $BIN/gc_cov_annotate.pl \
# --blasttaxid rand.nt.1e-5.megablast \
# --assembly $ASSEMBLY --bam *.bam --out bowtie_blobplot.txt \
# --taxdump $TAXDUMP --taxlist genus family phylum





