# P2: extract node-splitting comparison-level results from the netsplit object
# Writes results/P2_nodesplit_gi.csv and results/P2_nodesplit_panc.csv for figures.
suppressMessages({library(data.table); library(meta); library(netmeta)})
# Resolve repository root (parent of scripts/ when launched via Rscript)
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
ROOT <- if (!is.na(.file) && nzchar(.file)) {
  normalizePath(file.path(dirname(.file), ".."))
} else if (dir.exists("results")) {
  normalizePath(".")
} else if (dir.exists(file.path("..", "results"))) {
  normalizePath("..")
} else {
  stop("Run from the repository root, or: Rscript scripts/<script>.R")
}
results_dir <- file.path(ROOT, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
base <- results_dir
DL <- fread(file.path(results_dir, "P2_clean_drug_level_data.csv"))
cat("netmeta version:", as.character(packageVersion("netmeta")), "\n")

fit <- function(col) {
  d <- copy(DL); d[, ev := get(col)]
  d[, ev_adj := ev + 0.5 * (ev == 0)]
  pw <- pairwise(treat = treatment, event = ev_adj, n = n, studlab = study,
                 data = as.data.frame(d), sm = "RR")
  netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
          studlab = pw$studlab, sm = "RR", random = TRUE,
          reference.group = "Placebo", details.chkmultiarm = FALSE)
}

for (nm in c("responders_gi", "responders_panc")) {
  cat("\n##########", nm, "##########\n")
  net <- fit(nm)
  ns <- netsplit(net)
  cr <- as.data.frame(ns$compare.random, stringsAsFactors = FALSE)
  if (!"comparison" %in% names(cr) && !is.null(rownames(cr))) {
    cr$comparison <- rownames(cr)
  }
  cat("comparison labels (all):\n"); print(cr$comparison)
  cat("\n-- compare.random (testable only) --\n")
  testable <- cr[!is.na(cr$p), , drop = FALSE]
  print(testable)
  cat("\n-- p-values extracted --\n")
  if (nrow(testable) > 0) {
    pv <- as.numeric(testable$p)
    names(pv) <- testable$comparison
    print(round(pv, 4))
    cat("min p:", format.pval(min(pv, na.rm = TRUE), digits = 4),
        " comparisons p<0.05:", sum(pv < 0.05, na.rm = TRUE), "of", sum(!is.na(pv)), "\n")
  } else {
    cat("WARNING: zero testable comparisons\n")
  }
  outfile <- if (identical(nm, "responders_gi")) {
    "P2_nodesplit_gi.csv"
  } else {
    "P2_nodesplit_panc.csv"
  }
  fwrite(as.data.table(testable), file.path(results_dir, outfile))
  cat("wrote", outfile, "rows:", nrow(testable), "\n")
}
cat("\nDONE\n")