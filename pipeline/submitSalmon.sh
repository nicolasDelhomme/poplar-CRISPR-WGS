#!/bin/bash -l
set -eu -o pipefail

# load helpers
# shellcheck disable=SC1091
source "${SLURM_SUBMIT_DIR:-$(pwd)}/functions.sh"

# variables
ACCOUNT=u2018015
CPUs=20
DO=0
FORCE=0
KEY=

# keep these for when we start working on 1.1
REF=$(realpath ../reference)
SAMPLE_NUMBER=21

# all possible combinations of indexes
declare -A INXs=(
    ["Potra2"]="indices/salmon/Potra02_transcripts_salmon-version-1-dot-10-dot-3"
    ["Potrs1"]="../../Populus-tremuloides/v1.1/indices/salmon/Potrs01-mRNA_salmon-version-1-dot-10-dot-3"
    ["Potrx1"]="../../Populus-tremula_X_Populus-tremuloides/v1.0.1/indices/salmon/Potrx01_v1.0.1_gene-wo-introns_salmon-v-1-dot-10-dot-1"
    ["Potrx2"]="../../Populus-tremula_X_Populus-tremuloides/v2.0/indices/salmon/genome.mRNA.w.putative.pseudogene_with-decoy_salmon-version-1-dot-10-dot-3"
)

# usage
# shellcheck disable=SC2034
USAGETXT=\
"
  $0 [options]

  Options:
  -c number of CPUs to use, defaults to ${CPUs}
  -d do, not just print the commands
  -f force aligning
  -k redo a specific set
  -r alternate reference path, defaults to ${REF}
"

# tool
singularity=$(realpath ../singularity/salmon-1.10.3.sif)

# setup
export SINGULARITY_BINDPATH="/mnt:/mnt"

# options
# process the arguments
## get the options
while getopts c:dfk:r: option
do
  case "$option" in
    c) CPUs=${OPTARG};;
    d) DO=1;;
    f) FORCE=1;;
    k) KEY=${OPTARG};;
    r) REF=${OPTARG};;
	\?) ## unknown flag
	    abort "Unknow option";;
  esac
done
shift $(("$OPTIND" - 1))

# base results dir
out=$(realpath "../data")

# fastq dir
fastq="${out}/RNASeq/preprocessed/sortmerna/fqout"

# job list
cmds=()

if [ -n "${KEY}" ]; then
  INXs=( [${KEY}]=${INXs[$KEY]} )
fi

# submit
for INX in "${!INXs[@]}"; do

  # skip if not defined
  if [ "${INXs[$INX]}" != "" ]; then

    # build the out dir
    outdir="${out}/${INX}"
    [[ ! -d $outdir ]] && mkdir -p "$outdir"

    # check the outdir
    nfiles=$(find "${outdir}/" -name "quant.sf" | wc -l)

    # run
    if [ ${FORCE} -eq 1 ] || [ "${nfiles}" -ne "${SAMPLE_NUMBER}" ]; then
      # shellcheck disable=SC2044
      for f in $(find "${fastq}" -type f -name "*_fwd.fq.gz"); do
        fnam=$(basename "${f/_fwd.fq.gz/}")
        cmds+=("sbatch -A ${ACCOUNT} -t 48:00:00 \
-n ${CPUs} -o ${outdir}/${fnam}.out -e ${outdir}/${fnam}.err \
runSalmon.sh -t ${CPUs} ${singularity} ${REF}/${INXs[$INX]} ${f} ${f/_fwd.fq.gz/_rev.fq.gz} ${outdir}
")
      done
    else
      echo "Results exist for ${INX}. Rerun with -k ${INX} -f if you want to rebuild it."
    fi
  fi
done

if [ $DO -eq 1 ]; then
    for i in $(seq 0 $(( ${#cmds[@]} -1 ))); do
        eval "${cmds[$i]}"
    done
else
    for i in $(seq 0 $(( ${#cmds[@]} -1 ))); do
      echo "${cmds[$i]}"
    done
fi
