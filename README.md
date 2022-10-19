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

```sh
pyclone-vi fit -i /Users/bwilk/Documents/Projects/PVP/clonal_population_detection/chorangioma-pyclonevi-input.tsv -o chorangioma-pyclonevi-fit.h5 -c 40 -d beta-binomial -r 10
pyclone-vi write-results-file -i chorangioma-pyclonevi-fit.h5 -o chorangioma-pyclonevi-clonal-prediction.tsv
```

Running of PyClone-VI after reading into the various parameters and a few reviews from the tool authors as well as the
review paper originally leading to use of PyClone-VI:

```sh
pyclone-vi fit -i /Users/bwilk/Documents/Projects/PVP/clonal_population_detection/chorangioma-pyclonevi-input.tsv -o chorangioma-pyclonevi-binomial-100restarts-fit.h5 -c 40 -d binomial -r 100
pyclone-vi write-results-file -i chorangioma-pyclonevi-binomial-100restarts-fit.h5 -o chorangioma-pyclonevi-binomial-100restarts-clonal-prediction.tsv
```
