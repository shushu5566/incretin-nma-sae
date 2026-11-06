# P2: corrected figures for v1.1
# Fig3A/B driven by results/P2_qdecomp_gi.csv and results/P2_nodesplit_gi.csv
# (produced by fix_p2_corrected_analysis.R); fallback regenerates from clean data.
suppressMessages({library(ggplot2); library(data.table); library(patchwork); library(meta); library(netmeta)})
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
out  <- file.path(results_dir, "figures"); dir.create(out, showWarnings = FALSE)

gi <- fread(file.path(results_dir, "P2_nma_gi_corrected.csv"))
pc <- fread(file.path(results_dir, "P2_nma_pancreatitis_corrected.csv"))
DL <- fread(file.path(results_dir, "P2_clean_drug_level_data.csv"))
n_trials <- uniqueN(DL$study)
n_arms <- nrow(DL)
n_agents <- uniqueN(DL$treatment)
th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face = "bold"))

# ---- Figure 2: serious GI forest -------------------------------------------
gi[, sig := lo > 1]
p2 <- ggplot(gi, aes(x = RR, y = reorder(treatment, RR))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.7) +
  geom_point(aes(colour = sig), size = 2.8) +
  scale_colour_manual(values = c("FALSE" = "grey30", "TRUE" = "#C62828"), guide = "none") +
  scale_x_log10() +
  labs(x = "Risk ratio vs placebo (log scale)", y = NULL,
       title = "Serious gastrointestinal adverse events",
       subtitle = sprintf("%d trials, %d arms, %d agents; tau^2 and I^2 from random-effects netmeta",
                          n_trials, n_arms, n_agents)) + th
ggsave(file.path(out, "Fig2_GI_forest.png"), p2, width = 7, height = 4.2, dpi = 300)

# ---- Figure 4: pancreatitis forest ----------------------------------------
pc[, sig := lo > 1]
panc_events <- sum(DL$responders_panc)
p4 <- ggplot(pc, aes(x = RR, y = reorder(treatment, RR))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.7) +
  geom_point(colour = "#1565C0", size = 2.8) +
  scale_x_log10() +
  labs(x = "Risk ratio vs placebo (log scale)", y = NULL,
       title = "Serious pancreatitis",
       subtitle = sprintf("%d events across %d trials; no interval excludes 1",
                          panc_events, n_trials)) + th
ggsave(file.path(out, "Fig4_pancreatitis_forest.png"), p4, width = 7, height = 4.2, dpi = 300)

# ---- helpers: load or regenerate consistency inputs -----------------------
ensure_consistency_inputs <- function() {
  ns_path <- file.path(results_dir, "P2_nodesplit_gi.csv")
  qd_path <- file.path(results_dir, "P2_qdecomp_gi.csv")
  need_ns <- !file.exists(ns_path)
  need_qd <- !file.exists(qd_path)
  if (!need_ns && !need_qd) return(invisible(TRUE))
  cat("Regenerating missing consistency inputs via netsplit/decomp.design...\n")
  d <- copy(DL)
  d[, ev := responders_gi]
  d[, ev_adj := ev + 0.5 * (ev == 0)]
  pw <- pairwise(treat = treatment, event = ev_adj, n = n, studlab = study,
                 data = as.data.frame(d), sm = "RR")
  net <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
                 studlab = pw$studlab, sm = "RR", random = TRUE,
                 reference.group = "Placebo", details.chkmultiarm = FALSE)
  if (need_ns) {
    ns <- netsplit(net)
    cr <- as.data.frame(ns$compare.random, stringsAsFactors = FALSE)
    testable <- cr[!is.na(cr$p), , drop = FALSE]
    fwrite(as.data.table(testable), ns_path)
    cat("wrote", ns_path, "rows:", nrow(testable), "\n")
  }
  if (need_qd) {
    dl <- decomp.design(net)
    qtab <- as.data.frame(dl$Q.decomp)
    qtab$component <- rownames(qtab)
    fwrite(as.data.table(qtab), qd_path)
    cat("wrote", qd_path, "\n")
  }
  invisible(TRUE)
}
ensure_consistency_inputs()

# ---- Figure 3: consistency panel ------------------------------------------
qd_raw <- fread(file.path(results_dir, "P2_qdecomp_gi.csv"))
# Expect component column; map Within/Between rows
pick_row <- function(pat) {
  i <- grep(pat, qd_raw$component, ignore.case = TRUE)[1]
  if (is.na(i)) stop("Q.decomp missing row matching: ", pat)
  qd_raw[i]
}
w <- pick_row("^Within")
b <- pick_row("^Between")
qd <- data.frame(
  component = c("Within designs", "Between designs"),
  Q = c(as.numeric(w$Q), as.numeric(b$Q)),
  df = c(as.numeric(w$df), as.numeric(b$df)),
  p = c(as.numeric(w$pval), as.numeric(b$pval)),
  stringsAsFactors = FALSE)
p3a <- ggplot(qd, aes(x = Q, y = component)) +
  geom_col(width = 0.45, fill = "#2E7D32") +
  geom_text(aes(label = sprintf("Q=%.2f (df=%d), P=%.3f", Q, df, p)), hjust = -0.05, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.55))) +
  labs(x = "Q statistic", y = NULL,
       title = "A. Design-based Q decomposition (serious GI adverse events)") + th

ns <- fread(file.path(results_dir, "P2_nodesplit_gi.csv"))
if (!"comparison" %in% names(ns) || !"p" %in% names(ns)) {
  stop("P2_nodesplit_gi.csv must contain columns comparison and p")
}
ns <- ns[!is.na(p)]
n_testable <- nrow(ns)
if (n_testable == 0) stop("P2_nodesplit_gi.csv has zero testable comparisons")
p3b <- ggplot(ns, aes(x = p, y = reorder(comparison, p))) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "#C62828") +
  geom_point(size = 2.6, colour = "#6A1B9A") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "P value for direct vs indirect difference", y = NULL,
       title = {
         n_sig <- sum(ns$p < 0.05, na.rm = TRUE)
         if (n_sig == 0)
           sprintf("B. Node-splitting (%d testable comparisons; none significant)", n_testable)
         else
           sprintf("B. Node-splitting (%d testable comparisons; %d with P < 0.05)", n_testable, n_sig)
       }) + th

p3 <- p3a / p3b + plot_annotation(tag_levels = NULL)
ggsave(file.path(out, "Fig3_consistency.png"), p3, width = 7, height = 5.4, dpi = 300)

# ---- Figure 1: study selection flow ---------------------------------------
flow <- data.frame(
  x = c(1, 1, 1, 1, 1), y = c(5, 4, 3, 2, 1),
  label = c(
    "173 completed trials identified (ClinicalTrials.gov API v2, 13 intervention terms)",
    "50 trials with posted adverse-event results    (123 had none)",
    "+ pre-specified panel of 16 landmark incretin outcome trials\n→ 65 unique trials after deduplication on NCT number",
    "30 trials excluded before NMA:\n4 whole-trial dropouts (all arms unmapped to an incretin/placebo node:\nNCT00976391, NCT02229227, NCT03940742, NCT05813912)\n26 lacked a placebo arm or an incretin arm after arm-level mapping\nPK / drug-interaction / hepatic-impairment and insulin–incretin arms removed at arm level\nrepeated within-subject timepoint arms collapsed",
    sprintf("%d trials → %d drug-level arms → %d agents", n_trials, n_arms, n_agents)))
p1 <- ggplot(flow, aes(x = x, y = y)) +
  geom_label(aes(label = label), hjust = 0, size = 3.1, label.size = 0.3,
             fill = c("#E8F5E9","#E8F5E9","#E3F2FD","#FFF3E0","#E8F5E9"),
             label.padding = unit(0.45, "lines")) +
  xlim(0.9, 3.0) + ylim(0.5, 5.6) +
  labs(title = "Figure 1. Study selection") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))
ggsave(file.path(out, "Fig1_study_selection.png"), p1, width = 7.2, height = 4.2, dpi = 300)

cat("figures written:\n"); print(list.files(out))
cat("Fig3B comparisons:\n"); print(ns[, .(comparison, p = round(p, 4))])
cat("DONE\n")