#!/bin/bash

set -eu

proj=u2023008
mail=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/alternative_splicing/MEME
infile=intron_consensus_acceptor_sequences.fa
out=$in/pictogram_logos/acceptor_introns_consensus
sif=$(realpath ../singularity/memesuite_5.5.2.sif)

if [ ! -d $out ]; then
	mkdir -p $out
fi

export APPTAINER_BINDPATH="/mnt:/mnt"

sbatch --mail-user=$mail -A $proj -o ../logs/meme_logo.out -e ../logs/meme_logo.err \
runMEMEMakeLogo.sh $sif $in/$infile $out
