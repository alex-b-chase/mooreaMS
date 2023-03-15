#!/bin/bash

cd /Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/gene-centric-approach/biosynthetic_potential

cd blat_output

cat total.KSdomains.faa total.Cdomains.faa > ../cdhit_OBU/totalbiosyn.faa

## remove all reads <50 aa
reformat.sh in=totalbiosyn.faa out=totalbiosyn.filtered.faa minlength=50 ignorejunk

## need the read names to stay long and use high clustering
cd-hit -i totalbiosyn.filtered.faa -o totalbiosyn -c 0.8 -d 0 -g 1

clstr2txt.pl totalbiosyn.clstr > totalbiosyn.clstr.txt

cp totalbiosyn totalbiosyn.repseq.faa

