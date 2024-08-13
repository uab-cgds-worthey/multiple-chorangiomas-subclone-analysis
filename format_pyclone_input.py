import argparse
import gzip
import os
from pathlib import Path
import sys
import pandas as pd
from pyranges import PyRanges


def is_valid_output_file(p, arg):
    if os.access(Path(os.path.expandvars(arg)).parent, os.W_OK):
        return os.path.expandvars(arg)
    else:
        p.error(f"Output file {arg} can't be accessed or is invalid!")


def is_valid_file(p, arg):
    if not Path(os.path.expandvars(arg)).is_file():
        p.error("The file '%s' does not exist!" % arg)
    else:
        return os.path.expandvars(arg)


def get_AD_index(format_str):
    fcols = format_str.split(":")
    for index, col in enumerate(fcols):
        if col == "AD":
            return index

    return None


def read_muts_from_vcf(vcf, sample_id, normal_id):
    mut_dict = dict()
    with open(vcf, "rt") if vcf.endswith(".vcf") else gzip.open(vcf, "rt") as snvs_fp:
        normal_sample_col = 9
        somatic_sample_col = 10
        for line in snvs_fp:
            if line.startswith("#"):
                if not line.startswith("#CHROM"):
                    continue
                else:
                    cols = line.rstrip().split("\t")
                    normal_sample_col = cols.index(normal_id)
                    somatic_sample_col = cols.index(sample_id)

            vcfcols = line.rstrip().split("\t")
            # skip multiallelic sites as CNVKit doesn't use them for calculation
            if "," in vcfcols[4]:
                continue

            ad_index = get_AD_index(vcfcols[8])
            if ad_index:
                # check if variant info pulled from somatic variant, else pull counts from the germline variant
                ad_info = vcfcols[somatic_sample_col].split(":")[ad_index].split(",")
                if ad_info[1] in [".", "0"]:
                    ad_info = vcfcols[normal_sample_col].split(":")[ad_index].split(",")
                var_tuple = (vcfcols[0], vcfcols[1], vcfcols[3], vcfcols[4])
                mut_dict[var_tuple] = {
                    "mutation_id": "_".join(var_tuple),
                    "ref_counts": int(ad_info[0]),
                    "alt_counts": int(ad_info[1]),
                    "Chromosome": vcfcols[0],
                    "Start": int(vcfcols[1]),
                    "End": int(vcfcols[1]),
                }

    return mut_dict


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="Merge CNV and somatic biallelic SNV information together for input into PyCloneVI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    PARSER.add_argument(
        "--input",
        help=(
            "TSV file listing Sample ID, Sex (eg. F or M), Somatic SNV VCF, CNVKit segmetrics.tsv filepath, and "
            "merged results output filepath (one sample per line)"
        ),
        type=lambda x: is_valid_file(PARSER, x),
        default="data/config/sample_config.tsv",
        metavar="\b",
    )

    PARSER.add_argument(
        "--output",
        help="output file path to write merged PyClone-VI input file to",
        default="data/pyclone-input/merged_snvs_cnvs_pyclone.txt",
        type=lambda x: is_valid_output_file(PARSER, x),
        metavar="\b",
    )

    ARGS = PARSER.parse_args(sys.argv[1:])

    # read in sample information
    sample_info = pd.read_csv(
        ARGS.input, delimiter="\t", index_col="sample_id", header=0
    )

    # convert CNV segments file as a PyRange for easy feature intersection lookup
    # convert bi-allelic SNVs info into PyRange as well for intersection with CNV segments
    print("Converting CNVs and SNVs to PyRanges...")
    sample_cnv_segments = dict()
    sample_snvs = dict()
    for sample_id, row in sample_info.iterrows():
        # process CNV segments
        tempdf = pd.read_csv(
            row["cnv_segments"], delimiter="\t", index_col=False, header=0
        )
        # rename columns to make pyranges object
        tempdf.rename(
            columns={"chromosome": "Chromosome", "start": "Start", "end": "End"},
            inplace=True,
        )
        sample_cnv_segments[sample_id] = PyRanges(tempdf)

        # process SNVs from VCF file
        sample_snvs[sample_id] = read_muts_from_vcf(
            row["somatic_vcf"], sample_id, row["normal_id"]
        )

    # create a list of unique variant entries to fill in SNV counts with zeros
    # i.e., PyClone-Vi requires that all samples have all variants listed even if the variant wasn't called in the
    # given sample so we need to add those variants to a sample with values of zero to prevent it from pre-filtering
    # them out and not using for analysis
    var_set = set()  # makes unique set of Chromosome, Start tuples
    for snvs_dict in sample_snvs.values():
        var_set.update(snvs_dict.keys())

    for variant_tuple in var_set:
        for snv_dict in sample_snvs.values():
            if variant_tuple not in snv_dict:
                snv_dict[variant_tuple] = {
                    "mutation_id": "_".join(variant_tuple),
                    "ref_counts": 0,
                    "alt_counts": 0,
                    "Chromosome": variant_tuple[0],
                    "Start": int(variant_tuple[1]),
                    "End": int(variant_tuple[1]),
                }

    # join the CNV and SNV PyRanges for each sample (i.e. match overlapping CNV segment with a SNV)
    # see https://github.com/Roth-Lab/pyclone-vi#input-format for description of the PyClone-VI format we are writing out
    print("Joining CNV and SNV PyRanges for each sample and formatting columns...")
    sample_overlaps = dict()
    sex_chroms = {"chrX", "X", "chrY", "Y"}
    mut_df = None
    for sample, variants in sample_snvs.items():
        # join the PyRanges, then just work with the Pandas dataframe
        joined_df = (
            PyRanges(pd.DataFrame.from_dict(variants.values()))
            .join(sample_cnv_segments[sample], suffix="_cnv")
            .df
        )
        # drop chrY values b/c they confound analysis when cohort has a mix of males and females
        # Get indexes where name column has value john
        indexNames = joined_df[
            (joined_df["Chromosome"] == "chrY") | (joined_df["Chromosome"] == "Y")
        ].index
        # Delete these row indexes from dataFrame
        joined_df.drop(indexNames, inplace=True)
        # add the sample ID as a column
        joined_df["sample_id"] = sample
        # add normal copy number values, 2 if Female, 1 on sex chromosomes if Male
        if sample_info.at[sample, "sex"] == "F":
            joined_df["normal_cn"] = 2
        else:
            joined_df["normal_cn"] = joined_df.apply(
                lambda row: 1 if row["Chromosome"] in sex_chroms else 2, axis=1
            )
        # rename columns needed for output
        joined_df.rename(columns={"cn1": "major_cn", "cn2": "minor_cn"}, inplace=True)
        # only keep rows where allele specific copy number info is available
        joined_df = joined_df[joined_df["major_cn"].notna()]
        # ensure copy numbers are repesented in the output as ints
        joined_df["major_cn"] = joined_df["major_cn"].astype(int)
        joined_df["minor_cn"] = joined_df["minor_cn"].astype(int)

        # append or setup results
        if mut_df is None:
            mut_df = joined_df
        else:
            mut_df = pd.concat([mut_df, joined_df], ignore_index=True)

    print("Writing output...")
    columns = [
        "mutation_id",
        "sample_id",
        "ref_counts",
        "alt_counts",
        "major_cn",
        "minor_cn",
        "normal_cn",
    ]
    mut_df.sort_values(by=["mutation_id", "sample_id"], inplace=True)
    mut_df.to_csv(ARGS.output, sep="\t", index=False, columns=columns)

    print("Done. Outputs have been formatted and merged!")
