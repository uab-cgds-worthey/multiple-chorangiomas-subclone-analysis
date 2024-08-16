#!/usr/bin/env bash
SCRIPT_DIR="$(realpath $(dirname -- "${BASH_SOURCE[0]:-$0}"))"
# INPUT=${SCRIPT_DIR}/data/pyclone-input/merged_snvs_cnvs_pyclone.txt
INPUT=${SCRIPT_DIR}/data/pyclone-input/40depth_20vaf_merged_snvs_cnvs_pyclone.txt
OUTDIR=${SCRIPT_DIR}/data/40depth-20vaf-pyclone-output
mkdir -p ${OUTDIR}

set -eo pipefail

echo "Running pyclone fit"
pyclone-vi fit -i ${INPUT} -o ${OUTDIR}/pvp-pyclonevi-betabi-fit.h5 -c 40 -d beta-binomial -r 100 &>${OUTDIR}/pyclone-fit-betabi.log &
pyclone-vi fit -i ${INPUT} -o ${OUTDIR}/pvp-pyclonevi-binomial-fit.h5 -c 40 -d binomial -r 100 &>${OUTDIR}/pyclone-fit-binomial.log &

wait

echo "Writing results from pyclone"
pyclone-vi write-results-file -i ${OUTDIR}/pvp-pyclonevi-betabi-fit.h5 -o ${OUTDIR}/pvp-pyclonevi-betabi-clonal-prediction.tsv
pyclone-vi write-results-file -i ${OUTDIR}/pvp-pyclonevi-binomial-fit.h5 -o ${OUTDIR}/pvp-pyclonevi-binomial-clonal-prediction.tsv
echo "Done"
