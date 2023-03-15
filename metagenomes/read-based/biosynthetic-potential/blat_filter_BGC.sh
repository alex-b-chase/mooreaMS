#!/bin/bash

# take a metagenomic file, parse it for marker genes in the database

BASEDIR=/oasis/tscc/scratch/abchase/moorea/assembly/singlereadPHYLO
OUTDIR=$BASEDIR/biosyntheticP

rm -rf $OUTDIR
mkdir -p $OUTDIR

# using a new method to reduce computation time, substituting BLAT for BLAST
# with this, don't need to partition the reads into smaller chunks
# BLAT should be able to parse this in a few hours for each library from benchmarking on local machine

cd $BASEDIR/singlereadAA

for f in *.filter.total.faa
do

	output=${f%.filter.total.faa}

	for BGC in KS C 
	do
		echo "#!/bin/bash
#PBS -q condo
#PBS -N ${output}.${BGC}.blat
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:29:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

MGDIR=$BASEDIR/singlereadAA
BLASTDB=/oasis/tscc/scratch/abchase/refDB/BGC2KSC/napdosV2_${BGC}domains.faa
OUTDIR=$OUTDIR

cd \${MGDIR}

module load blat/35

blat -prot -fastMap -minIdentity=20 -out=blast8 \\
\$BLASTDB ${f} temp.${BGC}.${output}.txt

# subset the giant output file for only relevant information (i.e. query sequence and protein match)
cut -f1-2 temp.${BGC}.${output}.txt | \\
awk 'BEGIN{FS=\"\t\"; OFS=\"\t\"} {gsub(/^[^_]*_/, \"\", \$2); print}' | \\
sort -u > \$OUTDIR/${output}.blat.${BGC}.txt

rm -f temp.${BGC}.${output}.txt

# now subet the marker gene reads from the MG library
cut -f1 \$OUTDIR/${output}.blat.${BGC}.txt | \\
sort -u | sed 's/[ ]*$//' > \$OUTDIR/${output}.blat${BGC}.temp.txt

# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
filterbyname.sh \\
in=${f} \\
out=\$OUTDIR/${output}.blat.${BGC}markers.faa \\
names=\$OUTDIR/${output}.blat${BGC}.temp.txt \\
ow=t include=t tossjunk 2>/dev/null

rm -f \$OUTDIR/${output}.blat${BGC}.temp.txt


		" > $OUTDIR/$output.$BGC.sh

		qsub $OUTDIR/$output.$BGC.sh
		sleep 2s
	done
done

