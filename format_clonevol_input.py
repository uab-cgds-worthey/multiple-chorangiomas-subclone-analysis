import argparse
import csv
import os
import sys
from pathlib import Path


def is_valid_file(p, arg):
    if not Path(os.path.expandvars(arg)).is_file():
        p.error("The file '%s' does not exist!" % arg)
    else:
        return os.path.expandvars(arg)


def is_valid_output_dir(p, arg):
    dirpath = Path(os.path.expandvars(arg))
    if os.access(dirpath, os.W_OK) and dirpath.is_dir():
        return os.path.expandvars(arg)
    else:
        p.error(f"Output directory {arg} can't be accessed or is invalid!")


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="Join PyClone-VI output, variant info, and gene info for easy use with ClonEvol",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    PARSER.add_argument(
        "--input",
        help=(
            "TSV file listing sample info and clonal model file paths (one sample per line)"
        ),
        type=lambda x: is_valid_file(PARSER, x),
        default="data/config/clonevol_config.tsv",
        metavar="\b",
    )

    PARSER.add_argument(
        "--neatvars",
        help=(
            "TSV file listing variants marked as interesting in from variant analysis"
        ),
        type=lambda x: is_valid_file(PARSER, x),
        default="data/config/vars_of_interest.tsv",
        metavar="\b",
    )

    PARSER.add_argument(
        "--output",
        help="output directory path to write merged PyClone-VI input files to",
        default="data/clonevol-input/",
        type=lambda x: is_valid_output_dir(PARSER, x),
        metavar="\b",
    )

    PARSER.add_argument(
        "--sort_sample",
        help="sample to sort and renumber samples by, based on CCF of clusters for samples",
        default="CHORANGIOMA",
        type=str,
        metavar="\b",
    )

    ARGS = PARSER.parse_args(sys.argv[1:])

    # read in sample information
    sample_info = []
    with open(ARGS.input, "rt") as cmfp:
        reader = csv.DictReader(cmfp, delimiter="\t")
        sample_info = list(reader)

    # iterate samples and binomial and beta-binomial clustering output to format input for clonevol
    for model in ["binomial_clonal", "beta_binomial_clonal"]:
        variant_data = {}
        first_sample = sample_info[0]

        # read in the cluster and cancer cellular fraction (CCF) information output by pyclonevi
        with open(first_sample[model], "rt") as cmfp:
            reader = csv.DictReader(cmfp, delimiter="\t")
            for row in reader:
                mut = row["mutation_id"]
                cluster = row["cluster_id"]
                smp = row["sample_id"]

                if mut in variant_data:
                    variant_data[mut][f"{smp}_ccf"] = row["cellular_prevalence"]
                else:
                    variant_data[mut] = {
                        "mutation_id": mut,
                        "cluster": int(cluster) + 1,
                        f"{smp}_ccf": row["cellular_prevalence"],
                        "is_interesting": "FALSE",
                        "gene": "",
                    }

        # add variant allele fraction (VAF) for variants for each sample
        with open(first_sample["pyclone_input_vars"], "rt") as cmfp:
            reader = csv.DictReader(cmfp, delimiter="\t")
            for row in reader:
                mut = row["mutation_id"]
                if mut not in variant_data:
                    continue

                smp = row["sample_id"]

                depth = int(row["ref_counts"]) + int(row["alt_counts"])
                alt_cnts = int(row["alt_counts"])
                variant_data[mut][f"{smp}_vaf"] = (
                    0 if alt_cnts == 0 else float(alt_cnts) / float(depth)
                )

        # add interesting variants flag and gene information where present in sample
        with open(ARGS.neatvars, "rt") as cmfp:
            reader = csv.DictReader(cmfp, delimiter="\t")
            for row in reader:
                mut_id = (
                    row["chromosome"]
                    + "_"
                    + row["position"]
                    + "_"
                    + row["ref"]
                    + "_"
                    + row["alt"]
                )
                if mut_id in variant_data:
                    variant_data[mut_id]["is_interesting"] = "TRUE"
                    variant_data[mut_id]["gene"] = row["Gene"]

        # An issue with the model generation function in ClonEvol is that it requires reordering the clusters if
        # the founding clone isn't the one with the ID = 1.
        # See https://github.com/hdng/clonevol/issues/7#issuecomment-382786214
        # to resolve this easily we can sort clusters by CCF and relabel them by decreasing CCF per cluster
        uniq_ccf = {
            variant[f"{ARGS.sort_sample}_ccf"]: float(
                variant[f"{ARGS.sort_sample}_ccf"]
            )
            for variant in variant_data.values()
        }
        new_clust_num = 1
        new_clust_labels = {}
        for ccf in sorted(uniq_ccf.items(), key=lambda item: item[1], reverse=True):
            new_clust_labels[ccf[0]] = new_clust_num
            new_clust_num += 1

        for variant in variant_data.values():
            variant["cluster"] = new_clust_labels[variant[f"{ARGS.sort_sample}_ccf"]]

        # generate output file for the given model
        output_file = Path(ARGS.output) / f"pvp-{model}-clonevol.tsv"
        with output_file.open("wt") as outfp:
            fieldnames = list(next(iter(variant_data.values())).keys())
            writer = csv.DictWriter(outfp, fieldnames, delimiter="\t")

            writer.writeheader()
            writer.writerows(variant_data.values())

    print("Done. Outputs have been formatted and merged!")
