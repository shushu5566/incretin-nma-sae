# P2: repair node-splitting and design-based Q decomposition reporting
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
  cat("\n########## ", nm, " ##########\n")
  net <- fit(nm)
  cat("k =", net$k, " n =", net$n, " m =", net$m, "\n")
  cat("Q =", round(net$Q, 2), " df =", net$df.Q, " p =", format.pval(net$pval.Q, digits=4),
      " tau2 =", round(net$tau^2, 4), " I2 =", paste0(round(100*net$I2,1), "%"), "\n")

  cat("\n--- design-based decomposition ---\n")
  dd <- tryCatch(decomp.design(net), error = function(e) NULL)
  if (!is.null(dd)) {
    cat("class:", class(dd), "\n")
    for (comp in c("Q.decomp","Q.within","Q.between","df.within","df.between","pval.between")) {
      v <- dd[[comp]]
      if (!is.null(v)) { cat("--", comp, "--\n"); print(v) }
    }
  } else cat("unavailable\n")

  cat("\n--- node splitting ---\n")
  ns <- tryCatch(netsplit(net), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ns)) {
    cat("class:", class(ns), " length:", length(ns), "\n")
    cat("names:", paste(names(ns), collapse=", "), "\n")
    df <- tryCatch(as.data.frame(ns), error = function(e) { cat("as.data.frame ERR\n"); NULL })
    if (!is.null(df)) { cat("nrow:", nrow(df), " cols:", paste(names(df), collapse=", "), "\n"); print(df, digits = 3) }
  }

  cat("\n--- network connectivity / loops ---\n")
  cat("number of independent loops (df of network):", net$df, "\n")
  cat("network is star-shaped (no loops) if df == 0\n")
}
cat("\nDONE\n")
