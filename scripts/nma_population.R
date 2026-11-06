# P2: population-stratified NMA (obesity/overweight vs T2D) of serious GI SAE
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
sink(file.path(DIR, "nma_population_log.txt"), split = FALSE)
on.exit({
  tryCatch({ while (sink.number() > 0) sink() }, error = function(e) invisible(NULL))
}, add = TRUE)

# ---- build COMBINED drug-level dataset (pilot + expanded, dedup) ----
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
pilot <- read.csv(file.path(DIR, "nma_gi_sae_data.csv"), stringsAsFactors = FALSE)
pilot$drug <- sapply(pilot$treatment, drug_of)
pilot2 <- pilot %>% filter(!is.na(drug)) %>% group_by(study, drug) %>%
  summarise(responders = sum(responders), sampleSize = sum(sampleSize), .groups = "drop")
expanded <- read.csv(file.path(DIR, "nma_gi_druglevel.csv"), stringsAsFactors = FALSE)
names(expanded)[2] <- "drug"
dat <- bind_rows(pilot2, expanded) %>% distinct(study, drug, .keep_all = TRUE)
dat$treatment <- dat$drug
cat("combined:", nrow(dat), "arms,", length(unique(dat$study)), "trials\n")

titles <- read.delim(file.path(DIR, "ctgov_trial_list.txt"), sep = "\t", header = FALSE,
                     col.names = c("nct", "title"), stringsAsFactors = FALSE)

classify <- function(t) {
  tl <- tolower(t)
  if (grepl("obes|overweight|body weight|weight loss|weight management|bmi", tl)) return("Obesity")
  if (grepl("diabet|glycemic|glycaemic", tl)) return("T2D")
  if (grepl("cardiovascular", tl)) return("T2D")  # CVOTs are T2D populations
  "Other"
}
# pilot labels -> population map
pilot_pop <- c(STEP1="Obesity", STEP2="Obesity", STEP3="Obesity", SELECT="Obesity",
               "SUSTAIN-6"="T2D", "PIONEER-6"="T2D", FLOW="T2D", LEADER="T2D",
               SCALE="Obesity", "SURPASS-2"="T2D", "SURMOUNT-1"="Obesity",
               "SURPASS-4"="T2D", "SURMOUNT-2"="Obesity", REWIND="T2D",
               EXSCEL="T2D", "AMPLITUDE-O"="T2D", "TRIUMPH-obesity"="Obesity",
               "TRIUMPH-T2D"="T2D", "ATTAIN-1"="Obesity", "REDEFINE-T2D"="T2D",
               "SURVO-obesity"="Obesity", "MAZDUTIDE-CN"="Obesity")

# the druglevel file uses study = NCT for expanded, label for pilot - map both
dat$pop <- NA
# pilot studies are labelled by trial label; but druglevel.csv study column for pilot rows came from pilot2 (study=label)
# reconstruct: match by title where possible
for (i in seq_len(nrow(dat))) {
  s <- dat$study[i]
  if (s %in% names(pilot_pop)) { dat$pop[i] <- pilot_pop[[s]]; next }
  t <- titles$title[titles$nct == s]
  if (length(t) == 1 && !is.na(t)) dat$pop[i] <- classify(t)
}
cat("population distribution:\n"); print(table(dat$pop, useNA = "ifany"))

for (grp in c("Obesity", "T2D")) {
  dd <- dat[dat$pop == grp, ]
  if (nrow(dd) == 0) next
  cat("\n=== ", grp, ": ", length(unique(dd$study)), " trials ===\n", sep = "")
  cc <- dd
  cc$responders <- cc$responders + 0.5 * (cc$responders == 0)
  # drop single-arm trials within group
  armcnt <- table(cc$study)
  cc <- cc[cc$study %in% names(armcnt[armcnt >= 2]), ]
  if (length(unique(cc$study)) < 2) { cat("insufficient\n"); next }
  pw <- pairwise(treat = treatment, event = responders, n = sampleSize,
                 studlab = study, data = cc, sm = "RR")
  net <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
                 studlab = pw$studlab, sm = "RR", random = TRUE,
                 reference.group = "Placebo", details.chkmultiarm = FALSE)
  sm <- summary(net)$random
  res <- data.frame(treatment = rownames(sm$TE),
                    RR = exp(sm$TE[, "Placebo"]), lo = exp(sm$lower[, "Placebo"]),
                    hi = exp(sm$upper[, "Placebo"]), p = sm$p[, "Placebo"])
  res <- res[res$treatment != "Placebo" & !is.na(res$RR), ][order(res$RR), ]
  print(res)
  write.csv(res, file.path(DIR, paste0("nma_", tolower(grp), "_results.csv")), row.names = FALSE)
  cat("inconsistency Q p:", net$pval.Q, "\n")
}
cat("DONE\n")
sink()
