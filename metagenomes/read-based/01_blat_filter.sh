#!/bin/bash
#PBS -q condo
#PBS -N split2prod
#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:59:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

THREAD=4
BASEDIR=/oasis/tscc/scratch/abchase/moorea
OUTDIR=$BASEDIR/assembly

cd $OUTDIR

for f in *.filter.total.fa
do

	sampleID=${f%.filter.total.fa}

	partition.sh \
	in=$f \
	out=$sampleID.%.filter.total.fna ways=20

	i=0

	while [ $i -le 19 ]
	do
		echo "#!/bin/bash
#PBS -q condo
#PBS -N ${sampleID}_${i}_prod
#PBS -l nodes=1:ppn=4
#PBS -l walltime=2:59:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

REF=$sampleID
BASEDIR=${BASEDIR}
OUTDIR=\$BASEDIR/assembly
BLASTDB=/oasis/tscc/scratch/abchase/refDB/commmarkers/cmblastDB/total_markers.faa

cd \$OUTDIR

prodigal -i \$REF.${i}.filter.total.fna \\
-a \$REF.${i}.filter.total.faa -q \\
-f gff -p meta > \$REF.${i}.gff

rm -f \$REF.${i}.gff
rm -f \$REF.${i}.filter.total.fna
mv \$REF.${i}.filter.total.faa \$OUTDIR/singlereadPHYLO

cd \$OUTDIR/singlereadPHYLO

module load blat/35

blat -prot -fastMap -minIdentity=20 -out=blast8 \\
\$BLASTDB \$REF.${i}.filter.total.faa temp.20.\$REF.${i}.txt

# subset the giant output file for only relevant information (i.e. query sequence and protein match)
cut -f1-2 temp.20.\$REF.${i}.txt | \\
awk 'BEGIN{FS=\"\t\"; OFS=\"\t\"} {gsub(/^[^_]*_/, \"\", \$2); print}' | \\
sort -u > \$REF.${i}.blat.total.txt

rm -f temp.20.\$REF.${i}.txt

# now subet the marker gene reads from the MG library
cut -f1 \$REF.${i}.blat.total.txt | sort -u | sed 's/[ ]*$//' > \$REF.${i}.blat.temp.txt

# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
filterbyname.sh \\
in=\$REF.${i}.filter.total.faa \\
out=\$REF.${i}.blat.markers.faa \\
names=\$REF.${i}.blat.temp.txt ow=t include=t tossjunk

rm -f \$REF.${i}.blat.temp.txt


		" > $OUTDIR/singlereadPHYLO/${sampleID}_${i}.sh 
		qsub $OUTDIR/singlereadPHYLO/${sampleID}_${i}.sh 

		i=$(( $i + 1 ))
	done 
done

