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

## Copy-Number-aware clonal clustering using PyClone-VI

[PyClone-VI](https://github.com/Roth-Lab/pyclone-vi) is a faster update of PyClone which is used to infer clonal
population structure from genome sequencing data. It uses genome-wide CNV information as well as biallelic somatic SNVs
(and varaint allele fractions) to infer the clonal architecure of the tissue sample.

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

### Run PyClone-VI

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

[ClonEvol](https://github.com/hdng/clonevol) is a tool frequently used to visualize and inspect the results of clonal
evolution predictions. It's written in R and was easiest to use within an R markdown doc for tighter integration of the
analysis with notes.

Notes and the associated ClonEvol analysis code are all contained in [clonal-evo.Rmd](clonal-evo.Rmd). Further
information on this part of the analysis is contained within that document.

Final analysis results are contained within [data/clonevol-output](data/clonevol-output) containing the clonal evolution
models predicted from the somatic variant data.

### R and library versions

All required libraries are defined at the top of the R markdown and can be run to setup the analysis but the specific
versions used for producing output in this repo are defined below:

```
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.1 (2022-06-23)
 os       macOS 14.2.1
 system   x86_64, darwin17.0
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Chicago
 date     2024-01-09
 rstudio  2023.03.0+386 Cherry Blossom (desktop)
 pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version date (UTC) lib source
 cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
 callr         3.7.3   2022-11-02 [1] CRAN (R 4.2.0)
 cli           3.6.1   2023-03-23 [1] CRAN (R 4.2.0)
 clonevol    * 0.99.11 2022-11-02 [1] Github (hdng/clonevol@7aff737)
 colorspace    2.0-3   2022-02-21 [1] CRAN (R 4.2.0)
 crayon        1.5.1   2022-03-26 [1] CRAN (R 4.2.0)
 devtools      2.4.5   2022-10-11 [1] CRAN (R 4.2.0)
 digest        0.6.33  2023-07-07 [1] CRAN (R 4.2.0)
 dplyr         1.1.2   2023-04-20 [1] CRAN (R 4.2.0)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
 evaluate      0.23    2023-11-01 [1] CRAN (R 4.2.0)
 fansi         1.0.4   2023-01-22 [1] CRAN (R 4.2.0)
 farver        2.1.1   2022-07-06 [1] CRAN (R 4.2.0)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
 fs            1.5.2   2021-12-08 [1] CRAN (R 4.2.0)
 generics      0.1.3   2022-07-05 [1] CRAN (R 4.2.0)
 ggplot2     * 3.4.4   2023-10-12 [1] CRAN (R 4.2.0)
 glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
 gridBase    * 0.4-7   2014-02-24 [1] CRAN (R 4.2.0)
 gridExtra   * 2.3     2017-09-09 [1] CRAN (R 4.2.0)
 gtable        0.3.1   2022-09-01 [1] CRAN (R 4.2.0)
 htmltools     0.5.3   2022-07-18 [1] CRAN (R 4.2.0)
 htmlwidgets   1.5.4   2021-09-08 [1] CRAN (R 4.2.0)
 httpuv        1.6.5   2022-01-05 [1] CRAN (R 4.2.0)
 igraph      * 1.5.1   2023-08-10 [1] CRAN (R 4.2.0)
 knitr         1.39    2022-04-26 [1] CRAN (R 4.2.0)
 labeling      0.4.2   2020-10-20 [1] CRAN (R 4.2.0)
 later         1.3.0   2021-08-18 [1] CRAN (R 4.2.0)
 lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.2.0)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
 memoise       2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
 mime          0.12    2021-09-28 [1] CRAN (R 4.2.0)
 miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.2.0)
 munsell       0.5.0   2018-06-12 [1] CRAN (R 4.2.0)
 packcircles * 0.3.6   2023-09-08 [1] CRAN (R 4.2.0)
 pillar        1.9.0   2023-03-22 [1] CRAN (R 4.2.0)
 pkgbuild      1.3.1   2021-12-20 [1] CRAN (R 4.2.0)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
 pkgload       1.3.3   2023-09-22 [1] CRAN (R 4.2.0)
 prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
 processx      3.8.3   2023-12-10 [1] CRAN (R 4.2.0)
 profvis       0.3.7   2020-11-02 [1] CRAN (R 4.2.0)
 promises      1.2.0.1 2021-02-11 [1] CRAN (R 4.2.0)
 ps            1.7.5   2023-04-18 [1] CRAN (R 4.2.0)
 purrr         1.0.1   2023-01-10 [1] CRAN (R 4.2.0)
 R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
 Rcpp          1.0.9   2022-07-08 [1] CRAN (R 4.2.0)
 remotes       2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
 rlang         1.1.1   2023-04-28 [1] CRAN (R 4.2.0)
 rmarkdown     2.14    2022-04-25 [1] CRAN (R 4.2.0)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.2.0)
 scales        1.2.1   2022-08-20 [1] CRAN (R 4.2.0)
 sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
 shiny         1.7.2   2022-07-19 [1] CRAN (R 4.2.0)
 stringi       1.7.12  2023-01-11 [1] CRAN (R 4.2.0)
 stringr       1.5.0   2022-12-02 [1] CRAN (R 4.2.0)
 tibble        3.2.1   2023-03-20 [1] CRAN (R 4.2.0)
 tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.2.0)
 urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.2.0)
 usethis       2.1.6   2022-05-25 [1] CRAN (R 4.2.0)
 utf8          1.2.3   2023-01-31 [1] CRAN (R 4.2.0)
 vctrs         0.6.2   2023-04-19 [1] CRAN (R 4.2.0)
 withr         2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
 xfun          0.32    2022-08-10 [1] CRAN (R 4.2.1)
 xtable        1.8-4   2019-04-21 [1] CRAN (R 4.2.0)
 yaml          2.3.5   2022-02-21 [1] CRAN (R 4.2.0)

 [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library

────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```