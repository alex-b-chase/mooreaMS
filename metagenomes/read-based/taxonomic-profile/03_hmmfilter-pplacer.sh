#!/bin/bash

### THIS WILL APPLY A SECONDARY FILTER FOR THE MG READS FOR CURTO CORE GENES
### THE INPUT WILL BE THE FILTERED BLATp READS ${protein}.blat.faa
### NEED TO RUN THE MARKERS_SUBSET.SH FIRST!!!

# so each gene will most likely need an unique e-value
# we do not want to lose too much information so apply the loosest value to pass to PPLACER

REFDIR=/oasis/tscc/scratch/abchase/moorea/assembly/singlereadPHYLO
HMMPROF=/oasis/tscc/scratch/abchase/refDB/commmarkers/hmmprofiles
REFDB=/oasis/tscc/scratch/abchase/refDB/commmarkers/refpkg
OUTDIR=$REFDIR/pplacer

THREAD=4

cd $OUTDIR

# remove old files
rm -f *.e*
rm -f *.o*
rm -f *.jplace
rm -f *.csv
rm -f *.xml
rm -f *.log

# remove old files from HMMer scan earlier (not needed right now)

echo -e "coregene\te-valueHMMer" > $REFDIR/hmmerfiltered.txt

for blastout in *.blat.faa 
do 

	protein=${blastout%.blat.faa}

	echo "#!/bin/bash
#PBS -N ${protein}_pplacer
#PBS -q home-sio
#PBS -l nodes=1:ppn=${THREAD}
#PBS -l walltime=23:59:00

source ~/.bashrc

module load intel/2018.1.163
module load openmpi_ib/3.1.4

conda activate py36

PPLACERDIR=/home/abchase/software/pplacer

REFDIR=${REFDIR}
HMMPROF=${HMMPROF}
REFDB=${REFDB}

OUTDIR=\$REFDIR/pplacer 

cd \$OUTDIR

protein=${protein}

for minID in {15,20,25,30}
do

	cd \$OUTDIR

	hmmsearch --tblout \${protein}.hmm.txt -E 1e-\${minID} \\
	--cpu ${THREAD} \$HMMPROF/\${protein}p.hmm \${protein}.blat.faa > \${protein}.log.txt

	cat \${protein}.hmm.txt | cut -f1 -d' ' | sort | uniq > \${protein}.temp.txt

	# subset new HMMer filtered reads
	filterbyname.sh \\
	in=\${protein}.blat.faa \\
	out=total_\${protein}.\${minID}.hmm.faa \\
	names=\${protein}.temp.txt ow=t include=t tossjunk 2>/dev/null

	rm -f \${protein}.temp.txt
	rm -f \${protein}.hmm.txt
	rm -f \${protein}.log.txt

	cat total_\${protein}.\${minID}.hmm.faa | tr -d '*' > \$OUTDIR/\${protein}.\${minID}.temp.hmm.faa

	rm -f \$OUTDIR/\${protein}.\${minID}.fa

	clustalo --profile1 \$REFDB/\${protein}p.refpkg/\${protein}p.good.final.aln \\
	-i \${protein}.\${minID}.temp.hmm.faa -o \$OUTDIR/\${protein}.\${minID}.fa

	rm \${protein}.\${minID}.temp.hmm.faa

	### now can test input for pplacer 
	\$PPLACERDIR/pplacer --pretend \\
	-c \$REFDB/\${protein}p.refpkg \\
	\${protein}.\${minID}.fa > \${protein}.\${minID}.pplacer.log

	### check if pplacer worked
	if grep -Fxq \"everything looks OK.\" \${protein}.\${minID}.pplacer.log
	then
		echo -e \"\${protein}\\t\${minID}\" >> \$REFDIR/hmmerfiltered.txt
		break
	else 
		rm -f \$REFDIR/total_\${protein}.\${minID}.hmm.faa
		rm -f \$OUTDIR/\${protein}.\${minID}.fa
		continue
	fi
done

cd \$OUTDIR

rm -f \${protein}.*.pplacer.log

### now can run pplacer 
\$PPLACERDIR/pplacer -c \$REFDB/\${protein}p.refpkg \\
\${protein}.\${minID}.fa \\
-p --keep-at-most 20 

\$PPLACERDIR/guppy to_csv --point-mass --pp \${protein}.\${minID}.jplace > \${protein}.\${minID}.csv
\$PPLACERDIR/guppy fat --node-numbers --point-mass --pp \${protein}.*.jplace

conda deactivate

	" > $OUTDIR/${protein}.pplacer.sh

	qsub $OUTDIR/${protein}.pplacer.sh
	sleep 5
done




