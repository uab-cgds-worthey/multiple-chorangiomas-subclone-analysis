#!/usr/bin/env bash
SCRIPT_DIR="$(realpath $(dirname -- "${BASH_SOURCE[0]:-$0}"))"

usage() {
    cat <<-USAGE_MESSAGE
    usage: $0
    
    Run the clonal analysis of the chorangimoa and normal villi samples end-to-end with sepcified settings.

    options:
        -d | --depth    Min depth for somatic variant site to be included in analysis (default = 35)
        -v | --vaf      Min VAF for somatic variant site to be included in analysis (default = 0.08)
        -o | --out      Path to directory where output will be generated (default = ${SCRIPT_DIR}/data/output)
        -c | --config   Path to config TSV file with sample input information and metadata (default = ${SCRIPT_DIR}/data/metadata/sample_config.tsv)
        -h | --help     print usage info
USAGE_MESSAGE
}

while [ "$1" != "" ]; do
    case $1 in
    -d | --depth)
        shift
        DEPTH=$1
        ;;
    -v | --vaf)
        shift
        VAF=$1
        ;;
    -o | --out)
        shift
        OUTDIR=$1
        ;;
    -c | --config)
        shift
        CONFIG=$1
        ;;
    -h | --help)
        usage
        exit
        ;;
    *)
        usage
        exit 1
        ;;
    esac
    shift
done

if [[ -z $DEPTH ]]; then
    DEPTH=35
fi

if [[ -z $VAF ]]; then
    VAF=0.08
fi

if [[ -z $CONFIG ]]; then
    CONFIG=${SCRIPT_DIR}/data/metadata/sample_config.tsv
fi

if [[ -z $OUTDIR ]]; then
    OUTDIR=${SCRIPT_DIR}/data/output
fi

# enable conda to activate in a script, assumes conda env already made
eval "$(conda shell.bash hook)"

# set up other analysis paths and settings
METADATA_DIR=${SCRIPT_DIR}/data/metadata
CLONEVOL_CONFIG=${METADATA_DIR}/clonevol_input_config.tsv
INTERSTING_VARS=${METADATA_DIR}/vars_of_interest.tsv
FILT_CRITERIA=depth${DEPTH}-vaf${VAF}
PYCLONE_IN_DIR=${OUTDIR}/pyclone-${FILT_CRITERIA}-input
PYCLONE_OUT_DIR=${OUTDIR}/pyclone-${FILT_CRITERIA}-output
CLONEVOL_IN_DIR=${OUTDIR}/clonevol-${FILT_CRITERIA}-input
CLONEVOL_OUT_DIR=${OUTDIR}/clonevol-${FILT_CRITERIA}-output
mkdir -p ${PYCLONE_IN_DIR} ${PYCLONE_OUT_DIR} ${CLONEVOL_IN_DIR} ${CLONEVOL_OUT_DIR}

set -eo pipefail

echo "Formatting sample variant information as PyClone-VI input"
conda activate merge4pyclone
python ${SCRIPT_DIR}/src/format_pyclone_input.py --input ${CONFIG} --output ${PYCLONE_IN_DIR} --depth ${DEPTH} --vaf ${VAF}
conda deactivate

# compose clonevol config info along with running pyclone-vi
echo -e "sample_id\tpyclone_input_vars\tbinomial_clonal\tbeta_binomial_clonal" >${CLONEVOL_CONFIG}

for PYCLONE_INPUT in ${PYCLONE_IN_DIR}/*_merged_variant_info.tsv; do
    smp=$(echo ${PYCLONE_INPUT} | xargs basename | cut -d'_' -f1)

    echo "Running pyclone fit for ${smp}"
    conda activate pyclone-vi

    bifitfile=${PYCLONE_OUT_DIR}/${smp}-pyclonevi-binomial-fit.h5
    betabifitfile=${PYCLONE_OUT_DIR}/${smp}-pyclonevi-betabi-fit.h5
    pyclone-vi fit -i ${PYCLONE_INPUT} -o ${betabifitfile} -c 40 -d beta-binomial -r 50 &>${PYCLONE_OUT_DIR}/pyclone-fit-betabi.log &
    pyclone-vi fit -i ${PYCLONE_INPUT} -o ${bifitfile} -c 40 -d binomial -r 50 &>${PYCLONE_OUT_DIR}/pyclone-fit-binomial.log &

    wait

    echo "Writing results from pyclone for ${smp}"
    PYCLONE_BETA_PRED=${PYCLONE_OUT_DIR}/${smp}-pyclonevi-betabi-clonal-prediction.tsv
    PYCLONE_BINOMIAL_PRED=${PYCLONE_OUT_DIR}/${smp}-pyclonevi-binomial-clonal-prediction.tsv
    pyclone-vi write-results-file -i ${bifitfile} -o ${PYCLONE_BINOMIAL_PRED}
    pyclone-vi write-results-file -i ${betabifitfile} -o ${PYCLONE_BETA_PRED}
    echo -e "${smp}\t${PYCLONE_INPUT}\t${PYCLONE_BINOMIAL_PRED}\t${PYCLONE_BETA_PRED}" >>${CLONEVOL_CONFIG}
    conda deactivate
    echo "Done"
    echo ""
done

echo "Formatting input for ClonEvol"
conda activate merge4pyclone
python ${SCRIPT_DIR}/src/format_clonevol_input.py --input ${CLONEVOL_CONFIG} --output ${CLONEVOL_IN_DIR} --neatvars ${INTERSTING_VARS}
conda deactivate

echo "Running ClonEvol"
Rscript --vanilla ${SCRIPT_DIR}/src/clonevol_model.R -i ${CLONEVOL_IN_DIR} -o ${CLONEVOL_OUT_DIR}

echo "Analysis complete!"
