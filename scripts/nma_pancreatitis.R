# P2: pancreatitis drug-level NMA on combined dataset
suppressMessages({library(dplyr); library(meta); library(netmeta)})
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
sink(file.path(DIR, "nma_pancreatitis_log.txt"), split = FALSE)
on.exit({
  tryCatch({ while (sink.number() > 0) sink() }, error = function(e) invisible(NULL))
}, add = TRUE)

# build combined dataset with pancreatitis outcome
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
# pilot pancreatitis (from ctgov_nma_data.csv - has pancreatitis_affected)
pilot_raw <- read.csv(file.path(DIR, "ctgov_nma_data.csv"), stringsAsFactors = FALSE)
pilot_raw$drug <- sapply(pilot_raw$arm, function(a) {
  d <- drug_of(gsub("Semaglutide", "Sema", a)); d
})
# fix: use the node mapping logic from build_nma_data
node_map <- c("Semaglutide 2.4: Treatment Period"="Sema", "Semaglutide 2.4 mg"="Sema",
              "Semaglutide 1.0 mg"="Sema1.0", "Semaglutide 0.5 mg"="Sema0.5",
              "Oral Semaglutide"="OralSema", "Semaglutide"="Sema",
              "Liraglutide 3.0 mg, Pre-diabetes"="Lira", "Liraglutide 3.0 mg, no Pre-diabetes"="Lira",
              "Liraglutide"="Lira", "5 mg Tirzepatide"="Tirz", "10 mg Tirzepatide"="Tirz",
              "15 mg Tirzepatide"="Tirz", "Tirzepatide (MTD)"="Tirz",
              "1.5 mg Dulaglutide"="Dula", "Dulaglutide"="Dula", "Exenatide Once Weekly"="Exe",
              "Efpeglenatide 4 mg"="Efpe", "Efpeglenatide 6 mg"="Efpe",
              "6 mg Orforglipron"="Orfor", "12 mg Orforglipron"="Orfor", "36 mg Orforglipron"="Orfor")
pilot_raw$node <- node_map[pilot_raw$arm]
pilot_raw$node[grepl("Placebo", pilot_raw$arm)] <- "Placebo"
pilot_p <- pilot_raw %>% filter(!is.na(node)) %>%
  mutate(drug = sapply(node, function(n) if (n == "Placebo") "Placebo" else drug_of(n))) %>%
  filter(!is.na(drug)) %>%
  group_by(study = label, drug) %>%
  summarise(responders = sum(pancreatitis_affected), sampleSize = sum(arm_n), .groups = "drop")

# expanded pancreatitis (ctgov_expanded_data.csv has pancreatitis_affected)
exp_raw <- read.csv(file.path(DIR, "ctgov_expanded_data.csv"), stringsAsFactors = FALSE)
# reuse curate mapping by importing the same logic - simplified: rerun curation for pancreatitis
import_r <- function() {}
exp_raw$node <- sapply(exp_raw$arm, function(a) {
  x <- a
  if (grepl("placebo", x, ignore.case = TRUE)) return("Placebo")
  m <- regexpr("semaglutide", x, ignore.case = TRUE)
  if (m > 0) return("Sema")
  if (grepl("liraglutide", x, ignore.case = TRUE)) return("Lira")
  if (grepl("tirzepatide|ly3298176", x, ignore.case = TRUE)) return("Tirz")
  if (grepl("dulaglutide|ly2189265", x, ignore.case = TRUE)) return("Dula")
  if (grepl("exenatide", x, ignore.case = TRUE)) return("Exe")
  if (grepl("lixisenatide", x, ignore.case = TRUE)) return("Lixi")
  if (grepl("albiglutide", x, ignore.case = TRUE)) return("Albig")
  if (grepl("efpeglenatide", x, ignore.case = TRUE)) return("Efpe")
  if (grepl("ly3502970|orforglipron", x, ignore.case = TRUE)) return("Orfor")
  if (grepl("ly3437943|retatrutide", x, ignore.case = TRUE)) return("Reta")
  if (grepl("456906|survodutide", x, ignore.case = TRUE)) return("Survo")
  NA
})
exp_p <- exp_raw %>% filter(!is.na(node)) %>%
  mutate(drug = sapply(node, function(n) if (n == "Placebo") "Placebo" else drug_of(n))) %>%
  filter(!is.na(drug)) %>%
  group_by(study = nct, drug) %>%
  summarise(responders = sum(pancreatitis_affected), sampleSize = sum(arm_n), .groups = "drop")

combined <- bind_rows(pilot_p, exp_p) %>% distinct(study, drug, .keep_all = TRUE)
cat("combined pancreatitis: arms:", nrow(combined), "| trials:", length(unique(combined$study)), "\n")

cc <- combined
cc$responders <- cc$responders + 0.5 * (cc$responders == 0)
pw <- pairwise(treat = drug, event = responders, n = sampleSize, studlab = study, data = cc, sm = "RR")
net <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
               studlab = pw$studlab, sm = "RR", random = TRUE, reference.group = "Placebo",
               details.chkmultiarm = FALSE)
sm <- summary(net)$random
res <- data.frame(treatment = rownames(sm$TE),
                  RR = exp(sm$TE[, "Placebo"]), lo = exp(sm$lower[, "Placebo"]),
                  hi = exp(sm$upper[, "Placebo"]), p = sm$p[, "Placebo"])
res <- res[res$treatment != "Placebo", ][order(res$RR), ]
print(res)
write.csv(res, file.path(DIR, "nma_pancreatitis_results.csv"), row.names = FALSE)
cat("inconsistency Q p:", net$pval.Q, "\n")
cat("DONE\n")
sink()
