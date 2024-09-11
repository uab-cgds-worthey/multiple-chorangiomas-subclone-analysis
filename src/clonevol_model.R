# A simple script to efficiently run ClonEvol on multiple samples
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory containing clonal assignments, VAFs, and driver designations" # nolint
)

parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to deposit models into"
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if (is.null(args$input_dir) || is.null(args$output_dir)) {
  write("Input and Output directory paths must be specified!\n", stderr()) # nolint
  stop("Incorrect Input")
}

library(clonevol)

if (!dir.exists(args$output_dir)) dir.create(args$output_dir)

cluster_files <- list.files(args$input_dir, pattern = ".+\\.tsv", full.names = TRUE) # nolint
nonint_cols <- c("mutation_id", "cluster", "is_interesting", "gene")

# define clone colors
clone_colors <- NULL

for (cluster_tsv in cluster_files) {
  base_tsv_name <- basename(cluster_tsv)
  cat("Generating models for ", base_tsv_name, "\n")
  # determin colors for different models
  if (grepl("CHORANGIOMA", base_tsv_name) && grepl("beta", base_tsv_name)) {
    clone_colors <- c("#648FFF", "#009E73")
  } else {
    clone_colors <- c("#648FFF", "#CC79A7", "#009E73")
  }

  df <- read.table(file = cluster_tsv, header = TRUE, sep = "\t")

  # if only 1 cluster identified then skip b/c there is no model to be made
  if (length(unique(df$cluster)) < 2) {
    cat("Only a single clone detected for ", base_tsv_name, "\n")
    next
  }

  df[, !colnames(df) %in% nonint_cols] <- df[, !colnames(df) %in% nonint_cols] * 100 # nolint

  # ensure ordering of clusters is accending numerically to use with ClonEvol
  df <- df[order(df$cluster), ]

  # create dummy sample b/c ClonEvol no longer works with a single sample
  df$D.vaf <- df$vaf

  # rename VAF columns for better sample labeling in output diagrams
  smp_id <- unlist(strsplit(base_tsv_name, "-"))[1]
  vaf_col_names <- grep("vaf", colnames(df), value = TRUE)
  sample_names <- gsub("vaf", smp_id, vaf_col_names)
  sample_names <- gsub("NORM", "Normal_Villi", sample_names)
  sample_names <- gsub("CHORANGIOMA", "Chorangioma", sample_names)
  df[, sample_names] <- df[, vaf_col_names]
  vaf_col_names <- sample_names

  # compute clonal models
  y <- infer.clonal.models(
    variants = df,
    cluster.col.name = "cluster",
    vaf.col.names = vaf_col_names,
    clone.colors = clone_colors,
    cancer.initiation.model = "monoclonal",
    subclonal.test = "bootstrap",
    subclonal.test.model = "non-parametric",
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = "mean",
    ignore.clusters = NULL,
    min.cluster.vaf = 0.01,
    # min probability that CCF(clone) is non-negative
    sum.p = 0.05,
    # alpha level in confidence interval estimate for CCF(clone)
    alpha = 0.05
  )

  # map interesting variants/genes onto the tree if present
  if (dim(df[df$is_interesting, ])[1] > 0) {
    y <- transfer.events.to.consensus.trees(
      y,
      df[df$is_interesting, ],
      cluster.col.name = "cluster",
      event.col.name = "gene"
    )
  }

  # scale size of # of variants in cluster by square root transform
  y <- convert.consensus.tree.clone.to.branch(y, branch.scale = "log2")

  # plot the models
  plot.clonal.models(
    y,
    # box plot parameters
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = "grey50",
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = " Variant Allele Fraction",
    fancy.variant.boxplot.sample.title.size = 10,
    # bell plot parameters
    clone.shape = "bell",
    bell.event = TRUE,
    bell.event.label.color = "blue",
    bell.event.label.angle = 60,
    clone.time.step.scale = 1.25,
    bell.curve.step = 1,
    # node-based consensus tree parameters
    merged.tree.plot = FALSE,
    # branch-based consensus tree parameters
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ",",
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 0.25,
    mtcab.node.size = 4,
    mtcab.node.label.size = 1.5,
    mtcab.node.text.size = 1.4,
    # cellular population parameters
    cell.plot = TRUE,
    num.cells = 200,
    cell.border.size = 0.25,
    cell.border.color = "black",
    clone.grouping = "horizontal",
    # meta-parameters
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    # output figure parameters
    out.dir = args$output_dir,
    out.format = "pdf",
    out.prefix = sub("_clonal-clonevol.tsv", "-clonal-model", base_tsv_name), # nolint
    overwrite.output = TRUE,
    width = 11,
    height = 6,
    # vector of width scales for each panel from left to right
    panel.widths = c(3, 3, 3, 2)
  )
}
