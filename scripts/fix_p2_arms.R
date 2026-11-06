# P2 step 2: enumerate arm labels to build a robust drug mapping
suppressMessages(library(data.table))
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
exp <- fread(file.path(base, "ctgov_expanded_data.csv"))
pil <- fread(file.path(base, "ctgov_nma_data.csv"))
cat("=== all distinct arm labels, EXPANDED (", length(unique(exp$arm)), ") ===\n")
print(sort(unique(exp$arm)))
cat("\n=== all distinct arm labels, PILOT (", length(unique(pil$arm)), ") ===\n")
print(sort(unique(pil$arm)))
cat("\n=== trials with >1 timepoint-style arms ===\n")
print(exp[grepl("Baseline Through|Weeks Through|Treatment Period|Extension", arm), .N, by = nct][N > 4])
cat("DONE\n")
