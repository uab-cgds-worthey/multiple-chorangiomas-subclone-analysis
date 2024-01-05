# PyClone-VI

Installed following directions in https://github.com/Roth-Lab/pyclone-vi

## Installation

### Requirements

- Git v2.0+
- Mamba or Anaconda3

### Retrieve pipeline source code

This repository use git submodules, which needs to be pulled when cloning. Go to the directory of your choice and run
the command below.

```sh
git clone -b master git@github.com:uab-cgds-worthey/cnvkit_pipeline.git
```

Note that downloading this repository from GitHub, instead of cloning, may not fetch the submodules included.

### Create conda environment

The conda environment will install all necessary dependencies to run the various parts of this analysis.

```sh
# all commands assume your current working directory is the root directory of this repo

# create conda environments. Needed only the first time.
mamba env create -f env/input-formatter-env.yml
mamba env create -f env/pyclonevi-env.yml

# if you need to update the existing environment
mamba env update --file configs/envs/cnvkit.yaml
```

## Formatting inputs

Required Inputs:

- Mutect2 somatic bi-allelic SNVs with filter marked as `PASS`
- CNV segments file output from CNVKit called with the SNVs VCF listed above

Example sample config file for specifying merging criteria: [sample_config.tsv](.test/sample_config.tsv)

### Formmating SNVs

SNVs for clonal analysis were prefiltered for removal of FFPE artifact variants, `PASS` filter status from Mutect2, and
those where confidence in the call could be accertained (i.e., coverage and read support exceeded sequencing error
rate). Bialleleic SNV sites in the samples were obtained by using BCFTools following information in the
[BCFTools view command filter docs](https://samtools.github.io/bcftools/bcftools.html#view). Specific commands:

```sh
bcftools view -m2 -M2 -v snps -o LW001298.biallelic.mutect2.ideafix.PASS.hi-conf.vcf.gz LW001298.mutect2.ideafix.PASS.hi-conf.vcf.gz
bcftools view -m2 -M2 -v snps -o LW001299.biallelic.mutect2.ideafix.PASS.hi-conf.vcf.gz LW001299.mutect2.ideafix.PASS.hi-conf.vcf.gz
```

Resulting biallelic SNV VCF files are included in the [PyclonVI input dirctroy](data/pyclone-input).

### Merging inputs for PyCloneVI

To run PyClone-VI the somatic SNVs need to be combined with the CNV region information into the
[input format required by PyClone-VI](https://github.com/Roth-Lab/pyclone-vi/blob/master/README.md#input-format). The
script [format_pyclone_input.py](format_pyclone_input.py) script when run with the default parameters will process the
dataset for this analysis into the input format for PyClone-VI to run.

The formmating script can be run with the following

```sh
conda activate merge4pyclone
python format_pyclone_input.py
```

It'll produce the output file [merged_snvs_cnvs_pyclone.txt](data/pyclone-input/merged_snvs_cnvs_pyclone.txt) ready for
PyClone-VI analysis.

## Run PyClone-VI

Review article information indicated that the tool should be run with both the binomial and beta-binomial distributions
used for analysis. The number of clusters and restarts was recommended by a combination of review articles testing
multiple tools and documenation from the authors.

```sh
conda activate pyclone-vi
source run_pyclone.sh
```

PyClone-VI analysis will take several mintues to run and output the results in
[data/pyclone-output](data/pyclone-output). The \*.tsv files are the desired cluster predictions from both model types.
Based on information from the literature and the PyClone authors the binomial distrubution model is suitable for the
chorangioma tissue especially since we didn't detect any widespread CNV events.

## Visualizing Results with ClonEvol

[ClonEvol](https://github.com/hdng/clonevol) is a tool frequently used to visulalize and inspect the results of clonal
evolution predictions. It's written in R and was easiest to use within an R markdown doc for tighter integration of the
analysis with notes.

Notes and the associated ClonEvol analysis code are all contained in [clonal-evo.Rmd](clonal-evo.Rmd). Further
information on this part of the analysis is contained within that document.

Final analysis results are contained within [data/clonevol-output](data/clonevol-output) containing the clonal evolution
models predicted from the somatic variant data.

### R and library versions

All required libraries are defined at the top of the R markdown and can be run to setup the analysis but the specific
versions used for producing output in this repo are defined below:

- R version 4.2.1
  - attached base packages
    - stats graphics grDevices utils datasets methods base
  - other attached packages
    - clonevol_0.99.11 devtools_2.4.5 usethis_2.1.6
  - loaded via a namespace (and not attached)
    - tidyselect_1.2.0 xfun_0.32 remotes_2.4.2 purrr_1.0.1 generics_0.1.3 colorspace_2.0-3 vctrs_0.6.2 miniUI_0.1.1.1
    - htmltools_0.5.3 yaml_2.3.5 utf8_1.2.3 rlang_1.1.1 pkgbuild_1.3.1 later_1.3.0 urlchecker_1.0.1 pillar_1.9.0
    - glue_1.6.2 sessioninfo_1.2.2 lifecycle_1.0.3 stringr_1.5.0 munsell_0.5.0 gtable_0.3.1 htmlwidgets_1.5.4
      memoise_2.0.1
    - evaluate_0.23 knitr_1.39 callr_3.7.3 fastmap_1.1.0 httpuv_1.6.5 ps_1.7.5 curl_4.3.2 fansi_1.0.4
    - Rcpp_1.0.9 xtable_1.8-4 promises_1.2.0.1 scales_1.2.1 cachem_1.0.6 pkgload_1.3.3 mime_0.12 fs_1.5.2
    - ggplot2_3.4.4 digest_0.6.33 stringi_1.7.12 processx_3.8.3 dplyr_1.1.2 shiny_1.7.2 grid_4.2.1 cli_3.6.1
    - tools_4.2.1 magrittr_2.0.3 tibble_3.2.1 profvis_0.3.7 crayon_1.5.1 pkgconfig_2.0.3 ellipsis_0.3.2
      prettyunits_1.1.1
    - rmarkdown_2.14 rstudioapi_0.13 R6_2.5.1 compiler_4.2.1
