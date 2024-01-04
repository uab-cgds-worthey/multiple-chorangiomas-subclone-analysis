# PyClone-VI

Installed following directions in https://github.com/Roth-Lab/pyclone-vi

## Installation

### Requirements

-   Git v2.0+
-   Mamba or Anaconda3

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
conda env update --file configs/envs/cnvkit.yaml
```

## Formatting inputs

Required Inputs:

-   Mutect2 somatic bi-allelic SNVs with filter marked as `PASS`
-   CNV segments file output from CNVKit called with the SNVs VCF listed above

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
script [format_pyclone_input.py](format_pyclone_input.py) script when run with the default parameters will process
the dataset for this analysis into the input format for PyClone-VI to run. 

The formmating script can be run like

```sh
conda activate merge4pyclone
python format_pyclone_input.py
```

It'll produce the output file [merged_snvs_cnvs_pyclone.txt](data/pyclone-input/merged_snvs_cnvs_pyclone.txt) ready
for PyClone-VI analysis.

## Run PyClone-VI

Review article information indicated that the tool should be run with both the binomial and beta-binomial distributions
used for analysis. The number of clusters and restarts was recommended by a combination of review articles testing
multiple tools and documenation from the authors.

```sh
conda activate pyclone-vi
source run_pyclone.sh
```

PyClone-VI analysis will take several mintues to run the complete analyses and output the results
