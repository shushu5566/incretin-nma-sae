# P2 step 1: diagnose duplication in the source datasets
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
cat("=== ctgov_expanded_data.csv ===\n")
cat("rows:", nrow(exp), " unique nct:", length(unique(exp$nct)),
    " unique (nct,arm):", nrow(unique(exp[, .(nct, arm)])) , "\n")
cat("duplicated (nct,arm) rows:", sum(duplicated(exp[, .(nct, arm)])), "\n")
dups <- exp[duplicated(exp[, .(nct, arm)]) | duplicated(exp[, .(nct, arm)], fromLast = TRUE)]
if (nrow(dups)) { cat("\n--- duplicated (nct,arm) ---\n"); print(dups[order(nct, arm)]) }
cat("\nNCT01064687 rows:\n"); print(exp[nct == "NCT01064687"])
cat("\narms containing 'Baseline' or 'Week':\n")
print(head(exp[grepl("Baseline|Week|week|Month", arm), .(nct, arm, arm_n, pancreatitis_affected, serious_gi_affected)], 25))

pil <- fread(file.path(base, "ctgov_nma_data.csv"))
cat("\n=== ctgov_nma_data.csv (pilot) ===\n")
cat("rows:", nrow(pil), " unique nct:", length(unique(pil$nct)),
    " unique label:", length(unique(pil$label)), "\n")
cat("unique (nct,arm):", nrow(unique(pil[, .(nct, arm)])), "\n")
cat("labels:\n"); print(sort(unique(pil$label)))

cat("\n=== overlap between pilot and expanded, by NCT ===\n")
ov <- intersect(unique(pil$nct), unique(exp$nct))
cat("overlapping NCTs:", length(ov), "\n"); print(ov)
cat("\npilot NCTs NOT in expanded:", setdiff(unique(pil$nct), unique(exp$nct)), "\n")

cat("\n=== pancreatitis event totals ===\n")
cat("expanded: trials with any pancreatitis event:",
    length(unique(exp[pancreatitis_affected > 0, nct])), "of", length(unique(exp$nct)), "\n")
cat("expanded: arm rows with any pancreatitis event:", sum(exp$pancreatitis_affected > 0), "of", nrow(exp), "\n")
cat("expanded total pancreatitis events:", sum(exp$pancreatitis_affected), "\n")
cat("expanded total serious GI events:", sum(exp$serious_gi_affected), "\n")
cat("DONE\n")
