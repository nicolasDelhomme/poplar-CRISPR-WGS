#!/bin/bash -l

#safeguard

set -eu

# define variables
proj=u2023008
email=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/picard/v2.2
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/gatk/known_sites
ref=$(realpath ../reference/fasta/Potra02_genome_hardmasked.fasta)
gatk_sif=$(realpath ../singularity/gatk_4.2.6.1.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
for BAM in $(find $in -name "*.sorted_mkdup.bam"); do

sbatch -A $proj --mail-user $email -J 'gatk.'$(basename ${BAM/.sorted.bam/}) \
        runGenerateHQKnownSitesBQSR.sh $gatk_sif $BAM $ref $out
done
