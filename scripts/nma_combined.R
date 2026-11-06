# P2: combined drug-level NMA (flagship pilot + expanded) - serious GI SAE
suppressMessages({library(dplyr); library(meta); library(netmeta); library(gemtc)})
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
DIR <- results_dir
sink(file.path(DIR, "nma_combined_log.txt"), split = FALSE)
on.exit({
  tryCatch({ while (sink.number() > 0) sink() }, error = function(e) invisible(NULL))
}, add = TRUE)

# --- pilot dose-level -> drug-level ---
pilot <- read.csv(file.path(DIR, "nma_gi_sae_data.csv"), stringsAsFactors = FALSE)
drug_of <- function(node) {
  if (node == "Placebo") return("Placebo")
  for (prefix in c("Sema", "Lira", "Tirz", "Dula", "Exe", "Lixi", "Albig", "Efpe", "Orfor", "Reta", "Survo", "OralSema")) {
    if (startsWith(node, prefix)) {
      nm <- c(Sema="Semaglutide", Lira="Liraglutide", Tirz="Tirzepatide", Dula="Dulaglutide",
              Exe="Exenatide", Lixi="Lixisenatide", Albig="Albiglutide", Efpe="Efpeglenatide",
              Orfor="Orforglipron", Reta="Retatrutide", Survo="Survodutide", OralSema="OralSemaglutide")[prefix]
      return(nm)
    }
  }
  NA
}
pilot$drug <- sapply(pilot$treatment, drug_of)
pilot2 <- pilot %>% filter(!is.na(drug)) %>% group_by(study, drug) %>%
  summarise(responders = sum(responders), sampleSize = sum(sampleSize), .groups = "drop")

# --- expanded drug-level ---
expanded <- read.csv(file.path(DIR, "nma_gi_druglevel.csv"), stringsAsFactors = FALSE)
names(expanded)[2] <- "drug"

# combine, dedup by (study, drug) - pilot preferred (verified IDs)
combined <- bind_rows(pilot2, expanded) %>%
  distinct(study, drug, .keep_all = TRUE)

cat("combined: arms:", nrow(combined), "| trials:", length(unique(combined$study)),
    "| drugs:", length(unique(combined$drug)), "\n")
print(table(combined$drug))

# continuity correction for zero-event arms (pairwise RR)
cc <- combined
cc$responders <- cc$responders + 0.5 * (cc$responders == 0)
pw <- pairwise(treat = drug, event = responders, n = sampleSize, studlab = study, data = cc, sm = "RR")
net <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
               studlab = pw$studlab, sm = "RR", random = TRUE, reference.group = "Placebo",
               details.chkmultiarm = FALSE)
cat("\n=== Combined frequentist NMA (random, RR vs Placebo) ===\n")
sm <- summary(net)$random
res <- data.frame(treatment = rownames(sm$TE),
                  RR = exp(sm$TE[, "Placebo"]),
                  lo = exp(sm$lower[, "Placebo"]), hi = exp(sm$upper[, "Placebo"]),
                  p = sm$p[, "Placebo"])
res <- res[res$treatment != "Placebo", ]
res <- res[order(res$RR), ]
print(res)
write.csv(res, file.path(DIR, "nma_combined_results.csv"), row.names = FALSE)
cat("\ninconsistency Q p:", net$pval.Q, "\n")
png(file.path(DIR, "Fig_NMA_GI_combined.png"), width = 1600, height = 1200, res = 150)
forest(net, ref = "Placebo", digits = 2, main = "Serious GI SAE: combined frequentist NMA")
dev.off()
cat("DONE\n")
sink()
