#!/bin/bash

THREAD=4
BASEDIR=/oasis/tscc/scratch/abchase/moorea
RAWDIR=$BASEDIR/rawdata

cd $RAWDIR

# echo -e "sampleID\ttotalreads\tfilteredreads" > $BASEDIR/totalstats.txt

for f in M*_R1_001.fastq.gz
do

	sampleID=$(echo $f | cut -f1-2 -d'_')
	fread=$f
	rread=$(ls *${sampleID}*_R2_*)
	newname=$(echo ${sampleID} | sed s'/MO18_/MO_/g')

	# submit the HPC job as array job to BLASTp in smaller chunks
	echo "#!/bin/bash
#PBS -q condo
#PBS -N ${newname}_QC
#PBS -l nodes=1:ppn=${THREAD}
#PBS -l walltime=2:59:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

BASEDIR=${BASEDIR}
RAWDATADIR=\$BASEDIR/rawdata
OUTDIR=\$BASEDIR/assembly
BBMAPDIR=/home/abchase/software/bbmap/resources

REF=${newname}
THREAD=${THREAD}
FFILE=\$RAWDATADIR/${fread}
RFILE=\$RAWDATADIR/${rread}

cd \$OUTDIR

bbduk.sh qtrim=rl trimq=10 threads=${THREAD} \\
minlen=25 ktrim=r k=25 ref=\$BBMAPDIR/adapters.fa hdist=1 \\
in1=\$FFILE in2=\$RFILE \\
out1=\$REF.clean1.fq out2=\$REF.clean2.fq \\
tbo tpe &> \$REF.stats.txt

inputreads=\$(grep \"Input:\" \$REF.stats.txt | cut -f2 | cut -f1 -d' ' )
filteredreads=\$(grep \"Result:\" \$REF.stats.txt | cut -f2 | cut -f1 -d' ' )
rm -f \$REF.stats.txt

echo -e \"\${REF}\\t\${inputreads}\\t\${filteredreads}\" >> \$BASEDIR/totalstats.txt

repair.sh in=\$REF.clean1.fq in2=\$REF.clean2.fq \\
out=\$REF.filter.clean.R1.fq.gz out2=\$REF.filter.clean.R2.fq.gz overwrite=t

echo \"done!\"
sleep 10s

rm -f \$REF.clean1.fq
rm -f \$REF.clean2.fq

sleep 10s

	" > $BASEDIR/$newname.bbmap.sh

	# qsub $BASEDIR/$newname.bbmap.sh
	# sleep 15

done






