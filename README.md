# PyClone-VI

Installed following directions in https://github.com/Roth-Lab/pyclone-vi

## Formatting input

Required Inputs:

-   Mutect2 somatic bi-allelic SNVs with filter marked as `PASS`
-   CNV segments file output from CNVKit called with the SNVs VCF listed above

Example sample config file for specifying merging criteria: [sample_config.tsv](.test/sample_config.tsv)

### Formmating SNVs

SNVs for clonal analysis were prefiltered for removal of FFPE artifact variants, PASSing filter status from Mutect2, and
those where confidence in the call could be accertained (i.e., coverage and read support exceeded sequencing error
rate). Bialleleic SNV sites in the samples were obtained by using BCFTools following information in the
[BCFTools view command filter docs](https://samtools.github.io/bcftools/bcftools.html#view). Specific commands:

```sh
bcftools view -m2 -M2 -v snps -o LW001298.biallelic.mutect2.ideafix.PASS.hi-conf.vcf.gz LW001298.mutect2.ideafix.PASS.hi-conf.vcf.gz
bcftools view -m2 -M2 -v snps -o LW001299.biallelic.mutect2.ideafix.PASS.hi-conf.vcf.gz LW001299.mutect2.ideafix.PASS.hi-conf.vcf.gz
```

Resulting biallelic SNV VCF files are included in the [PyclonVI input dirctroy](data/pyclone-input).

### Merging inputs for PyCloneVI

```sh
conda activate merge4pyclone
```



## Running the tool

Tested with the following command to first explore if and how the tools worked and the output:

Directory where work with PyClone-VI is being done: /Users/bwilk/Documents/Projects/PVP/clonal_population_detection
**Data from analysis subject to move to permanent location eventually and this note will be updated**

Review article information indicated that the tool should be run with both the binomial and beta-binomial distributions
used for analysis. I'm still not sure if CNVs matter for Chorangioma analysis as there does not appear to be major CNV
alterations detected in the tissue. The number of clusters and restarts was recommended by a combination of review
articles testing multiple tools and documenation from the authors.

```sh
pyclone-vi fit -i chorangioma-pyclonevi-input.tsv -o chorangioma-pyclonevi-betabi-fit.h5 -c 10 -d beta-binomial -r 100
pyclone-vi fit -i chorangioma-pyclonevi-input.tsv -o chorangioma-pyclonevi-binomial-fit.h5 -c 10 -d binomial -r 100
pyclone-vi write-results-file -i chorangioma-pyclonevi-betabi-fit.h5 -o chorangioma-pyclonevi-betabi-clonal-prediction.tsv
pyclone-vi write-results-file -i chorangioma-pyclonevi-binomial-fit.h5 -o chorangioma-pyclonevi-binomial-clonal-prediction.tsv
```
