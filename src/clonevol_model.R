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

# Monkey patch the cross.rule.score() function in the clonevol package so that
# an error in scoring concensus models doesn't crash the application
# monkey patch for R came from this blog post https://dlukes.github.io/monkey-patching-in-r.html
# this resolves an issue resolved in this SO post https://stackoverflow.com/a/42036069/2892199
# `clone.ccf.combined.p` just needs to be initianlized before being set for some odd reson


#########################################################################
cross.rule.score <- function(x, meta.p.method = "fisher", exhaustive.mode = FALSE,
                             rank = TRUE, boot = NULL) {
  if (!is.null(x$matched) && x$num.matched.models > 0 && ncol(x$matched$index) > 1) {
    if (!is.null(x$matched$trimmed.merged.trees)) {
      cat("WARN: pruned trees found. pruneConsensusTrees must be rerun.\n")
    }
    samples <- names(x$models)
    num.models <- nrow(x$matched$index)
    x$matched$scores$max.clone.ccf.combined.p <- NA
    x$matched$clone.ccf.pvalues <- list()
    # foreach matched model, recalc score by combining p values across
    # samples for each clone
    for (i in 1:num.models) {
      trees <- NULL
      t <- x$match$merged.trees[[i]]
      p <- NULL
      for (s in samples) {
        mi <- x$models[[s]][[x$match$index[i, s]]]
        mi <- mi[!mi$excluded & !is.na(mi$parent), c("lab", "p.value")]
        colnames(mi) <- c("lab", s)
        if (is.null(p)) {
          p <- mi
        } else {
          p <- merge(p, mi, all = TRUE)
        }
      }
      # ppp <<- p
      if (ncol(p) == 2) { # single sample
        p$cmb.p <- apply(p[, c(2, 2)], 1, clonevol:::combine.p, method = meta.p.method)
      } else {
        p$cmb.p <- apply(p[, -1], 1, clonevol:::combine.p, method = meta.p.method)
      }
      # model score = max (combined p of each clone)
      x$matched$scores$max.clone.ccf.combined.p[i] <- max(p$cmb.p)
      # this is never used elsewhere and tanks the whole plotting
      x$matched$merged.trees[[i]]$clone.ccf.combined.p <- NA
      # save the whole pvalue matrix
      x$matched$clone.ccf.pvalues[[i]] <- p
    }
    # order matched models by new score
    idx <- seq(1, nrow(x$matched$scores))
    if (rank) {
      idx <- order(x$matched$scores$max.clone.ccf.combined.p)
    }
    x$matched$index <- x$matched$index[idx, , drop = FALSE]
    x$matched$scores <- x$matched$scores[idx, , drop = FALSE]
    x$matched$probs <- x$matched$probs[idx, , drop = FALSE]
    # order merged trees
    tmp <- list()
    for (i in idx) {
      tmp <- c(tmp, list(x$matched$merged.trees[[i]]))
    }
    x$matched$merged.trees <- tmp
    # order merged traces
    tmp <- list()
    for (i in idx) {
      tmp <- c(tmp, list(x$matched$merged.traces[[i]]))
    }
    x$matched$merged.traces <- tmp
    # order pvalues
    tmp <- list()
    for (i in idx) {
      tmp <- c(tmp, list(x$matched$clone.ccf.pvalues[[i]]))
    }
    x$matched$clone.ccf.pvalues <- tmp

    # remove previous obsolete model scores (which was very small probability)
    # x$matched$scores$model.prob = x$matched$scores$model.score
    # x$matched$scores$model.score = NULL
  }
  return(x)
}

assignInNamespace("cross.rule.score", cross.rule.score, "clonevol")
#########################################################################


if (!dir.exists(args$output_dir)) dir.create(args$output_dir)

cluster_files <- list.files(args$input_dir, pattern = ".+\\.tsv", full.names = TRUE) # nolint
nonint_cols <- c("mutation_id", "cluster", "is_interesting", "gene")

for (cluster_tsv in cluster_files) {
  base_tsv_name <- basename(cluster_tsv)
  cat("Generating models for ", base_tsv_name, "\n")

  df <- read.table(file = cluster_tsv, header = TRUE, sep = "\t")

  # if only 1 cluster identified then skip b/c there is no model to be made
  if (length(unique(df$cluster)) < 2) {
    cat("Only a single clone detected for ", base_tsv_name, "\n")
    next
  }

  df[, !colnames(df) %in% nonint_cols] <- df[, !colnames(df) %in% nonint_cols] * 100 # nolint

  # ensure ordering of clusters is accending numerically to use with ClonEvol
  df <- df[order(df$cluster), ]

  # rename VAF columns for better sample labeling in output diagrams
  vaf_col_names <- grep("_vaf", colnames(df), value = TRUE)
  sample_names <- gsub("_vaf", "", vaf_col_names)
  sample_names <- gsub("NORM_VILLI", "Normal_Villi", sample_names)
  sample_names <- gsub("CHORANGIOMA", "Chorangioma", sample_names)
  df[, sample_names] <- df[, vaf_col_names]
  vaf_col_names <- sample_names

  # define clone colors
  clone_colors <- NULL

  # compute clonal models
  y <- infer.clonal.models(
    variants = df,
    cluster.col.name = "cluster",
    vaf.col.names = vaf_col_names,
    clone.colors = clone_colors,
    cancer.initiation.model = "polyclonal",
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
