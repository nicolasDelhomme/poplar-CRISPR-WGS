#!/bin/bash -l

set -eu

# define variables
proj=u2023021
email=kristina.benevides@umu.se
in=/mnt/ada/projects/arabidopsis/swenkel/arabidopsis-frosty/bwa
out=/mnt/ada/projects/arabidopsis/swenkel/arabidopsis-frosty/picard
reference=$(realpath ../reference)
picard_sif=$(realpath ../singularity/picard_2.27.1.sif)
samtools_sif=$(realpath ../singularity/samtools_1.16.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
for BAM in $(find $in -name "*.sorted.bam"); do

sbatch -A $proj --mail-user $email -J $(basename ${BAM/.sorted.bam/}) \
	runPicardMarkDuplicatesWithMateCigar.sh $picard_sif $samtools_sif \
	$BAM $out
done
