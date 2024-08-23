#!/usr/bin/env bash
SCRIPT_DIR="$(realpath $(dirname -- "${BASH_SOURCE[0]:-$0}"))"

usage() {
    cat <<-USAGE_MESSAGE
    usage: $0
    
    Run the clonal analysis of the chorangimoa and normal villi samples end-to-end with sepcified settings.

    options:
        -d | --depth    Min depth for somatic variant site to be included in analysis (default = 35)
        -v | --vaf      Min VAF for somatic variant site to be included in analysis (default = 0.08)
        -i | --inc      Force inclusion of interesting somatic variants despite depth reqs (default = false)
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
    -i | --inc)
        shift
        INCVARS=1
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

if [[ -z $INCVARS ]]; then
    INCVARS=0
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
INTERSTING_VARS=${SCRIPT_DIR}/data/metadata/vars_of_interest.tsv
FILT_CRITERIA=depth${DEPTH}-vaf${VAF}
PYCLONE_IN_DIR=${OUTDIR}/pyclone-${FILT_CRITERIA}-input
PYCLONE_OUT_DIR=${OUTDIR}/pyclone-${FILT_CRITERIA}-output
CLONEVOL_IN_DIR=${OUTDIR}/clonevol-${FILT_CRITERIA}-input
CLONEVOL_OUT_DIR=${OUTDIR}/clonevol-${FILT_CRITERIA}-output
mkdir -p ${PYCLONE_IN_DIR} ${PYCLONE_OUT_DIR} ${CLONEVOL_IN_DIR} ${CLONEVOL_OUT_DIR}

set -eo pipefail

echo "Formatting sample variant information as PyClone-VI input"
PYCLONE_INPUT=${PYCLONE_IN_DIR}/merged_snvs_cnvs_pyclone.txt
conda activate merge4pyclone
if [ $INCVARS -eq 1 ]; then
    python ${SCRIPT_DIR}/src/format_pyclone_input.py --input ${CONFIG} --output ${PYCLONE_INPUT} --depth ${DEPTH} --vaf ${VAF} --incvars
else
    python ${SCRIPT_DIR}/src/format_pyclone_input.py --input ${CONFIG} --output ${PYCLONE_INPUT} --depth ${DEPTH} --vaf ${VAF}
fi

# echo "Running pyclone fit"
conda activate pyclone-vi
pyclone-vi fit -i ${PYCLONE_INPUT} -o ${PYCLONE_OUT_DIR}/pvp-pyclonevi-betabi-fit.h5 -c 40 -d beta-binomial -r 100 &>${PYCLONE_OUT_DIR}/pyclone-fit-betabi.log &
pyclone-vi fit -i ${PYCLONE_INPUT} -o ${PYCLONE_OUT_DIR}/pvp-pyclonevi-binomial-fit.h5 -c 40 -d binomial -r 100 &>${PYCLONE_OUT_DIR}/pyclone-fit-binomial.log &

wait

# echo "Writing results from pyclone"
PYCLONE_BETA_PRED=${PYCLONE_OUT_DIR}/pvp-pyclonevi-betabi-clonal-prediction.tsv
PYCLONE_BINOMIAL_PRED=${PYCLONE_OUT_DIR}/pvp-pyclonevi-binomial-clonal-prediction.txt
echo "Done"

# echo "Formatting input for ClonEvol"
conda activate merge4pyclone
CLONEVOL_CONFIG=${CLONEVOL_IN_DIR}/clonevol_input_config.tsv
echo -e "pyclone_input_vars\tbinomial_clonal\tbeta_binomial_clonal" >${CLONEVOL_CONFIG}
echo -e "${PYCLONE_INPUT}\t${PYCLONE_BINOMIAL_PRED}\t${PYCLONE_BETA_PRED}" >>${CLONEVOL_CONFIG}
python ${SCRIPT_DIR}/src/format_clonevol_input.py --input ${CLONEVOL_CONFIG} --output ${CLONEVOL_IN_DIR} --neatvars ${INTERSTING_VARS} --sort_sample CHORANGIOMA

echo "Running ClonEvol"
Rscript --vanilla ${SCRIPT_DIR}/src/clonevol_model.R -i ${CLONEVOL_IN_DIR} -o ${CLONEVOL_OUT_DIR}

echo "Analysis complete!"
