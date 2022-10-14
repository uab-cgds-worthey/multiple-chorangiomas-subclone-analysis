import argparse
import csv
import os
from pathlib import Path
import sys


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


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        description="Merge CNV gene information output from CNVkit for cbioportal loading",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    PARSER.add_argument(
        "--input",
        help="TSV file listing Sample ID, Somatic SNV VCF, and CNVKit segmetrics.tsv filepath to merge (one path per line)",
        type=lambda x: is_valid_file(PARSER, x),
        required=True,
        metavar="\b",
    )

    PARSER.add_argument(
        "--output",
        help="output file path to write merged PyClone-VI input file to",
        default="merged_snvs_cnvs_pyclone.txt",
        type=lambda x: is_valid_output_file(PARSER, x),
        metavar="\b",
    )

    ARGS = PARSER.parse_args(sys.argv[1:])

    # parse HGNC file for information to convert gene symbol to EntrezGene ID
    official_symbol_lookup = dict()
    previous_symbols_lookup = dict()
    with open(ARGS.hgnc) as fp:
        for index, line in enumerate(fp):
            if index > 0:
                cols = line.strip().split("\t")
                if len(cols) < 10:
                    continue
                official_symbol_lookup[cols[1]] = cols[9]
                for prev_symb in cols[4].split(","):
                    prev_symb = prev_symb.strip()
                    previous_symbols_lookup[prev_symb] = cols[9]

    # collect the unique list of genes listed in all files to be merged
    print("Collecting info...")
    genes = set()
    log2s = dict()
    with open(ARGS.cnvs) as fp:
        for line in fp:
            sample_info = line.split(",")
            sample_id = sample_info[0]
            log2s[sample_id] = dict()
            with open(sample_info[1].strip()) as log2fp:
                for dataline in log2fp:
                    if dataline.startswith("gene	chromosome"):
                        continue

                    cols = dataline.split("\t")
                    genes.add(cols[0])
                    log2s[sample_id][cols[0]] = cols[4]

    if ARGS.purecn:
        with open(ARGS.purecn) as fp:
            for line in fp:
                sample_info = line.split(",")
                sample_id = sample_info[0]
                log2s[sample_id] = dict()
                with open(sample_info[1].strip(), newline="") as log2fp:
                    reader = csv.DictReader(log2fp)
                    for row in reader:
                        genes.add(row["symbol"])
                        if row["cn_category"] == "Amplification":
                            log2s[sample_id][row["symbol"]] = 2
                        elif row["cn_category"] == "Gain":
                            log2s[sample_id][row["symbol"]] = 1
                        elif row["cn_category"] == "Neutral":
                            log2s[sample_id][row["symbol"]] = 0
                        elif row["cn_category"] == "Loss":
                            log2s[sample_id][row["symbol"]] = -1
                        else:
                            log2s[sample_id][row["symbol"]] = -2

    print("Merging continuous and descrete values...")
    # following guidelines listed in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8312221/ we can provide a descrete
    # value for CNV changes that may work better with cBioPortal via the following setup:
    # -> continuous values lower than −2.0 were listed as homozygous deletions (HOMDEL)
    # -> values between −2.0 and −1.0 were listed as hemizygous deletions (HETLOSS)
    # -> values between 1.0 and 2.0 were declared amplifications (GAIN),
    # -> values more than 2.0 multiple amplifications (AMP).
    # -> value was between −1.0 and 1.0 were declared diploid samples (DIPLOID), meaning they have no changes in the number of gene copies
    #
    # these values correspond with what cbioportal can use
    with open(ARGS.output_disc, "wt") as out_disc:
        with open(ARGS.output_cont, "wt") as out_cont:
            # print header information, cbioportal specific
            sample_keys = sorted(list(log2s.keys()))
            out_disc.write("Hugo_Symbol\tEntrez_Gene_Id\t" + "\t".join(sample_keys))
            out_cont.write("Hugo_Symbol\tEntrez_Gene_Id\t" + "\t".join(sample_keys))
            for gene in sorted(genes):
                out_cont.write("\n")
                out_disc.write("\n")
                gene_id = ""
                if gene in official_symbol_lookup.keys():
                    gene_id = official_symbol_lookup[gene]
                elif gene in previous_symbols_lookup.keys():
                    gene_id = previous_symbols_lookup[gene]

                out_disc.write(gene + "\t" + gene_id)
                out_cont.write(gene + "\t" + gene_id)
                for sample_id in sample_keys:
                    if gene in log2s[sample_id].keys():
                        out_cont.write("\t" + str(log2s[sample_id][gene]))
                        log2 = float(log2s[sample_id][gene])
                        if log2 <= -2:
                            out_disc.write("\t-2")
                        elif log2 > -2 and log2 <= -1:
                            out_disc.write("\t-1")
                        elif log2 > -1 and log2 < 1:
                            out_disc.write("\t0")
                        elif log2 >= 1 and log2 < 2:
                            out_disc.write("\t1")
                        else:
                            out_disc.write("\t2")

                    else:
                        out_cont.write("\tNA")
                        out_disc.write("\tNA")

    print(f"Done. Outputs have been merged!")
