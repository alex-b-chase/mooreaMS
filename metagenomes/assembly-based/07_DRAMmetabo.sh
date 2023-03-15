#!/bin/bash
#SBATCH --job-name=DRAM_MGs
#SBATCH -p medium-mem-1-m
#SBATCH --mem=512G
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=36    ## number of cores the job needs
#SBATCH --error=test-%J.err	## error log file
#SBATCH --output=test-%J.out	## output info file

source ~/.bashrc

conda activate DRAM

cd /work/users/abchase/moorea/assembly/binned_genomes/dRep

DRAM.py annotate \
-i 'dereplicated_genomes/*.fa' -o DRAMannotation \
--gtdb_taxonomy gtdb_class/gtdbtk.bac120.summary.tsv \
--checkm_quality metabat.checkm.tsv \
--threads 36

DRAM.py distill \
-i DRAMannotation/annotations.tsv -o DRAMdistill \
--trna_path DRAMannotation/trnas.tsv \
--rrna_path DRAMannotation/rrnas.tsv

conda deactivate
