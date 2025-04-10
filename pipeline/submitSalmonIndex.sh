#!/bin/bash -l
set -eu -o pipefail

# load helpers
# shellcheck disable=SC1091
source "${SLURM_SUBMIT_DIR:-$(pwd)}/functions.sh"

# variables
CPUs=20
DO=0
DMEM="--mem=240GB"
FORCE=0
KEY=
REF=$(realpath ../reference)

# all possible combinations of fasta
declare -A FASTAs=(
    ["Potra2"]="fasta/Potra02_transcripts.fasta.gz"
    ["Potrs1"]="../../Populus-tremuloides/v1.1/fasta/Potrs01-mRNA.fa.gz"
    ["Potrx1"]=""
    ["Potrx2"]=""
)

# usage
# shellcheck disable=SC2034
USAGETXT=\
"
  $0 [options]

  Options:
  -c number of CPUs to use, defaults to ${CPUs}
  -d do, not just print the commands
  -f force rebuilding
  -k rebuild a specific set
  -r alternate reference path, defaults to ${REF}
"

# tool
singularity=$(realpath ../singularity/salmon-1.10.3.sif)

# decoy
#decoy=$(realpath ../reference/fasta/Potra02_genome.fasta.gz)

# setup
export SINGULARITY_BINDPATH="/mnt:/mnt"

# sanity
# extend all with the ref and check
for FASTA in ${!FASTAs[@]}; do
    [[ "${FASTAs[${FASTA}]}" != "" ]] && [[ ! -f ${REF}/${FASTAs[${FASTA}]} ]] && abort "The file ${REF}/${FASTAs[${FASTA}]} does not exist"
    [[ "${FASTAs[${FASTA}]}" != "" ]] && FASTAs[${FASTA}]=${REF}/${FASTAs[${FASTA}]}
done

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

# job list
cmds=()

if [ -n "${KEY}" ]; then
  FASTAs=( [${KEY}]=${FASTAs[$KEY]} )
fi

# submit
for FASTA in "${!FASTAs[@]}"; do

  # skip if not defined
  if [ "${FASTAs[$FASTA]}" != "" ]; then

    # get the fasta
    fa="${FASTAs[$FASTA]}"
    fdir=$(basename "${fa/.f*a.gz//}")

    # build the out dir
    outdir="$(dirname "${FASTAs[${FASTA}]}")/../indices/salmon/${fdir}_salmon-version-1-dot-10-dot-3"
    [[ ! -d $outdir ]] && mkdir -p "$outdir"

    # run w/o decoy
    if [ ${FORCE} -eq 1 ] || [ ! -f "${outdir}/versionInfo.json" ]; then
      cmds+=("sbatch -A kaw -t 48:00:00 \
-n ${CPUs} -o ${outdir}.out -e ${outdir}.err \
runSalmonIndex.sh -t ${CPUs} ${singularity} ${fa} ${outdir}
")
    else
      echo "An index exists for ${FASTA}. Rerun with -k ${FASTA} -f if you want to rebuild it."
    fi

#     # rebuild the out dir
#     outdir="${out}/${fdir}_with-decoy_salmon-version-1-dot-10-dot-3"
#     [[ ! -d $outdir ]] && mkdir -p "$outdir"

#     # run with decoy
#     if [ ${FORCE} -eq 1 ] || [ ! -f "${outdir}/versionInfo.json" ]; then
#       cmds+=("sbatch -A kaw -t 48:00:00 ${DMEM} \
# -n ${CPUs} -o ${outdir}.out -e ${outdir}.err \
# runSalmonIndex.sh -d ${decoy} -t ${CPUs} ${singularity} ${fa} ${outdir}
# ")
    # else
    #   echo "An index with decoy exists for ${FASTA}. Rerun with -k ${FASTA} -f if you want to rebuild it."
    # fi
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
