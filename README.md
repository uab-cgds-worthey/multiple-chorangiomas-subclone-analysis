# PyClone-VI

Installed following directions in https://github.com/Roth-Lab/pyclone-vi

## Formatting input

Required Inputs:

-   Mutect2 somatic bi-allelic SNVs with filter marked as `PASS`
-   CNV segments file output from CNVKit called with the SNVs VCF listed above

Example sample config file for specifying merging criteria: [sample_config.tsv](.test/sample_config.tsv)

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
