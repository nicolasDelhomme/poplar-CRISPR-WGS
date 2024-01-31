#!/bin/bash
set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

proj=u2023008
genomedir=/mnt/ada/projects/aspseq/T89-assembly/data/
outdir=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/T89_STAR/
star=/mnt/picea/storage/singularity/kogia/star_2.7.10a.sif

mkdir -p $outdir/index
#zcat $genomedir/T89.hap1.fa.gz $genomedir/T89.hap2.fa.gz > $outdir/genome.fa
sbatch -A $proj ./runSTARGenomeGenerate.sh -n $star $outdir/index $outdir/genome.fa
