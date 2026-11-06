# =============================================================================
# P2 CORRECTED ANALYSIS
# Fixes:
#  (1) exclude PK / drug-interaction / hepatic-impairment studies (protocol excludes them)
#      -- the expanded dataset contained DDI arms (midazolam, simvastatin, ...) and
#         repeated within-subject orforglipron day-range arms
#  (2) collapse NCT01064687 timepoint arms to a single period (56 weeks)
#  (3) dedup pilot vs expanded on NCT (the pilot DOES carry nct; the old code keyed
#      the pilot by trial label, so dedup silently failed)
#  (4) exclude insulin-incretin combination arms (separate interventions)
#  (5) re-run netmeta with I^2, tau^2, Q and node-splitting on the FULL network
#  (6) report event sparsity and a SELECT-excluded sensitivity analysis
# =============================================================================
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
sink(file.path(results_dir, "P2_corrected_log.txt"))
on.exit({
  tryCatch({ while (sink.number() > 0) sink() }, error = function(e) invisible(NULL))
}, add = TRUE)

cat("netmeta version:", as.character(packageVersion("netmeta")), "\n")
cat("meta version:", as.character(packageVersion("meta")), "\n")

exp <- fread(file = file.path(base, "ctgov_expanded_data.csv"))
pil <- fread(file = file.path(base, "ctgov_nma_data.csv"))

# ---- titles (to identify PK/DDI studies) ------------------------------------
tit <- tryCatch(fread(file = file.path(base, "nma_gi_expanded.csv")), error = function(e) NULL)
if (!is.null(tit) && "title" %in% names(tit)) {
  tt <- unique(tit[, .(nct = study, title)])
  exp <- merge(exp, tt, by = "nct", all.x = TRUE)
} else exp[, title := NA_character_]

cat("=== STEP 1: trial-level exclusions ===\n")
pk_pat <- "pharmacokinetic|drug-drug|drug interaction|bioequivalence|hepatic impairment|relative bioavailability|food effect|QT|thorough"
title_map_present <- !is.null(tit) && "title" %in% names(tit) && any(!is.na(exp$title) & nzchar(exp$title))
if (!title_map_present) {
  cat("title-level PK/DDI exclusion SKIPPED (mapping table nma_gi_expanded.csv not in public package);\n")
  cat("arm-level regex exclusion (PK probe drugs / hepatic labels / insulin-incretin combos) still applied in STEP 4.\n")
  pk_trials <- character(0)
} else {
  pk_trials <- unique(exp[grepl(pk_pat, title, ignore.case = TRUE), nct])
  cat("PK / interaction / hepatic-impairment trials identified and EXCLUDED:", length(pk_trials), "\n")
  print(unique(exp[nct %in% pk_trials, .(nct, title)]))
}
exp2 <- exp[!(nct %in% pk_trials)]

# ---- STEP 2: collapse timepoint arms ---------------------------------------
cat("\n=== STEP 2: timepoint collapse (NCT01064687) ===\n")
cat("rows before:", nrow(exp2), "\n")
exp2 <- exp2[!(nct == "NCT01064687" & grepl("26 Weeks Through 56 Weeks", arm))]
exp2 <- exp2[!(nct == "NCT01064687" & grepl("Baseline Through 26 Weeks", arm))]
cat("rows after collapsing to the 56-week period:", nrow(exp2), "\n")

# ---- STEP 3: drug mapping ---------------------------------------------------
incretins <- c("semaglutide","liraglutide","tirzepatide","dulaglutide","exenatide",
               "lixisenatide","albiglutide","efpeglenatide","orforglipron",
               "retatrutide","survodutide","mazdutide")
pk_drugs <- "midazolam|acetaminophen|rosuvastatin|simvastatin|moxifloxacin|digoxin|sodium bicarbonate"

map_drug <- function(arm) {
  a <- tolower(arm)
  if (grepl(pk_drugs, a)) return(NA_character_)
  if (grepl("hepatic impairment|normal hepatic function", a)) return(NA_character_)
  hit <- NA_character_
  if (grepl("semaglutide", a))                        hit <- "Semaglutide"
  else if (grepl("liraglutide", a))                   hit <- "Liraglutide"
  else if (grepl("tirzepatide|ly3298176", a))         hit <- "Tirzepatide"
  else if (grepl("dulaglutide|ly2189265", a))         hit <- "Dulaglutide"
  else if (grepl("exenatide", a))                     hit <- "Exenatide"
  else if (grepl("lixisenatide|iglarlixi", a))        hit <- "Lixisenatide"
  else if (grepl("albiglutide", a))                   hit <- "Albiglutide"
  else if (grepl("efpeglenatide", a))                 hit <- "Efpeglenatide"
  else if (grepl("orforglipron|ly3502970", a))        hit <- "Orforglipron"
  else if (grepl("retatrutide|ly3437943", a))         hit <- "Retatrutide"
  else if (grepl("survodutide|bi ?456906", a))        hit <- "Survodutide"
  else if (grepl("mazdutide|ly2944876", a))           hit <- "Mazdutide"
  if (!is.na(hit) && grepl("insulin", a)) return(NA_character_)   # insulin-incretin combos
  if (is.na(hit) && grepl("placebo", a)) return("Placebo")
  hit
}

# ---- build the arm-level dataset (expanded + pilot), deduped on NCT ---------
mk <- function(d) {
  d <- as.data.table(d)
  if (!"title" %in% names(d)) d[, title := NA_character_]
  d[, drug := sapply(arm, map_drug)]
  d[, .(nct, arm, drug, arm_n, serious_gi_affected, pancreatitis_affected)]
}
A <- rbind(mk(exp2), mk(pil), fill = TRUE)
cat("\n=== STEP 4: dedup pilot vs expanded on NCT ===\n")
cat("arm rows before dedup:", nrow(A), " unique nct:", length(unique(A$nct)), "\n")
ov <- intersect(unique(pil$nct), unique(exp2$nct))
cat("NCTs present in both sources:", length(ov), paste(ov, collapse=", "), "\n")
A <- A[!is.na(nct)]
A <- unique(A, by = c("nct", "arm"))
cat("arm rows after dedup:", nrow(A), " unique nct:", length(unique(A$nct)), "\n")

n_unmap_rows <- sum(is.na(A$drug))
cat("\nunmapped arms (excluded by arm-level drug map):", n_unmap_rows, "rows\n")
print(sort(unique(A[is.na(drug), arm]))[1:25])

# Trials where EVERY arm fails mapping disappear entirely from the network
nct_before_map <- unique(A$nct)
whole_drop <- A[, .(all_na = all(is.na(drug))), by = nct][all_na == TRUE, nct]
cat("whole-trial dropouts (all arms unmapped):", length(whole_drop), "\n")
if (length(whole_drop)) print(whole_drop)

A <- A[!is.na(drug)]
# restrict to trials with a placebo arm and >=1 incretin arm
ok <- A[, .(has_pbo = any(drug == "Placebo"), has_drug = any(drug != "Placebo")), by = nct]
keep <- ok[has_pbo & has_drug, nct]
n_fail_design <- length(setdiff(unique(A$nct), keep))
cat("\ntrials with both a placebo and an incretin arm:", length(keep), "of", length(unique(A$nct)), "\n")
cat("excluded for no placebo or no incretin arm:", n_fail_design, "\n")
cat("flowchart exclusion total (whole-trial unmapped + design):",
    length(whole_drop) + n_fail_design, "of", length(nct_before_map), "\n")
A <- A[nct %in% keep]

# drug-level aggregation per (nct, drug)
DL <- A[, .(responders_gi = sum(serious_gi_affected),
            responders_panc = sum(pancreatitis_affected),
            n = sum(arm_n)), by = .(study = nct, treatment = drug)]
cat("\nfinal drug-level arms:", nrow(DL), " trials:", length(unique(DL$study)),
    " drugs:", length(unique(DL$treatment)), "\n")
cat("\narms per drug:\n"); print(DL[, .N, by = treatment][order(-N)])

cat("\n=== STEP 5: event sparsity ===\n")
cat("serious GI events total:", sum(DL$responders_gi),
    " in arms with >=1 event:", sum(DL$responders_gi > 0), "of", nrow(DL), "\n")
cat("pancreatitis events total:", sum(DL$responders_panc),
    " in arms with >=1 event:", sum(DL$responders_panc > 0), "of", nrow(DL), "\n")
cat("trials contributing >=1 pancreatitis event:",
    length(unique(DL[responders_panc > 0, study])), "of", length(unique(DL$study)), "\n")

# ---- helpers for netmeta 3.6.x API -----------------------------------------
extract_nodesplit <- function(ns) {
  # Prefer compare.random (netmeta 3.6.1); as.data.frame(netsplit) may fail.
  cr <- NULL
  if (!is.null(ns$compare.random)) {
    cr <- ns$compare.random
  } else {
    cr <- tryCatch(as.data.frame(ns), error = function(e) {
      cat("WARNING: as.data.frame(netsplit) failed:", conditionMessage(e), "\n")
      NULL
    })
  }
  if (is.null(cr)) return(NULL)
  cr <- as.data.frame(cr, stringsAsFactors = FALSE)
  if (!"comparison" %in% names(cr) && !is.null(rownames(cr))) {
    cr$comparison <- rownames(cr)
  }
  pcol <- intersect(c("p", "pval", "p.value"), names(cr))[1]
  if (is.na(pcol)) {
    cat("WARNING: no p-value column in netsplit result; columns:",
        paste(names(cr), collapse = ", "), "\n")
    return(NULL)
  }
  cr$p <- suppressWarnings(as.numeric(cr[[pcol]]))
  cr
}

extract_qdecomp <- function(dl) {
  # netmeta 3.6.1: results live in $Q.decomp (rows Total / Within / Between).
  # Older API used scalar $Q.within / $Q.between / $pval.between.
  out <- list(Q.within = NA_real_, df.within = NA_real_, p.within = NA_real_,
              Q.between = NA_real_, df.between = NA_real_, p.between = NA_real_,
              source = NA_character_)
  if (is.null(dl)) return(out)
  qd <- dl$Q.decomp
  if (!is.null(qd)) {
    qd <- as.data.frame(qd)
    rn <- rownames(qd)
    pick <- function(pat) {
      i <- grep(pat, rn, ignore.case = TRUE)[1]
      if (is.na(i)) return(NULL)
      list(Q = as.numeric(qd$Q[i]), df = as.numeric(qd$df[i]),
           p = as.numeric(qd$pval[i]))
    }
    w <- pick("^Within")
    b <- pick("^Between")
    if (!is.null(w)) {
      out$Q.within <- w$Q; out$df.within <- w$df; out$p.within <- w$p
    }
    if (!is.null(b)) {
      out$Q.between <- b$Q; out$df.between <- b$df; out$p.between <- b$p
    }
    out$source <- "Q.decomp"
    return(out)
  }
  if (!is.null(dl$Q.within) || !is.null(dl$Q.between)) {
    out$Q.within <- suppressWarnings(as.numeric(dl$Q.within)[1])
    out$Q.between <- suppressWarnings(as.numeric(dl$Q.between)[1])
    out$p.between <- suppressWarnings(as.numeric(dl$pval.between)[1])
    out$source <- "legacy_scalars"
  }
  out
}

# ---- NMA -------------------------------------------------------------------
run_nma <- function(outcome_col, label, drop_study = NULL, save_aux = FALSE) {
  d <- copy(DL)
  if (!is.null(drop_study)) d <- d[!(study %in% drop_study)]
  d[, ev := get(outcome_col)]
  d[, ev_adj := ev + 0.5 * (ev == 0)]     # documented continuity correction
  pw <- pairwise(treat = treatment, event = ev_adj, n = n, studlab = study,
                 data = as.data.frame(d), sm = "RR")
  net <- netmeta(TE = pw$TE, seTE = pw$seTE, treat1 = pw$treat1, treat2 = pw$treat2,
                 studlab = pw$studlab, sm = "RR", random = TRUE,
                 reference.group = "Placebo", details.chkmultiarm = FALSE)
  sm <- summary(net)$random
  res <- data.frame(treatment = rownames(sm$TE),
                    RR = exp(sm$TE[, "Placebo"]),
                    lo = exp(sm$lower[, "Placebo"]),
                    hi = exp(sm$upper[, "Placebo"]),
                    p  = sm$p[, "Placebo"])
  res <- res[res$treatment != "Placebo", ]
  res <- res[order(res$RR), ]
  cat("\n================ ", label, if (!is.null(drop_study)) paste0(" (excluding ", drop_study, ")") else "", " ================\n")
  # netmeta: k=studies, n=treatments, m=pairwise comparisons (NOT arm rows)
  cat("trials:", net$k, " treatments:", net$n, " pairwise comparisons:", net$m, "\n")
  cat("Q (total):", round(net$Q, 2), " df:", net$df.Q, " p:", format.pval(net$pval.Q, digits = 4), "\n")
  cat("tau^2:", round(net$tau^2, 4), "  tau:", round(net$tau, 4), "\n")
  cat("I^2 (total):", paste0(round(100 * net$I2, 1), "%"), "\n")
  print(res, row.names = FALSE, digits = 3)

  cat("\n-- node-splitting (full network) --\n")
  ns <- tryCatch(netsplit(net), error = function(e) {
    cat("ERROR: netsplit failed:", conditionMessage(e), "\n"); NULL
  })
  nsdf <- NULL
  if (!is.null(ns)) {
    nsdf <- extract_nodesplit(ns)
    if (is.null(nsdf)) {
      cat("WARNING: node-splitting object present but no extractable table\n")
    } else {
      testable <- nsdf[!is.na(nsdf$p), , drop = FALSE]
      cat("comparisons in table:", nrow(nsdf),
          " testable (non-NA p):", nrow(testable), "\n")
      if (nrow(testable) == 0) {
        cat("WARNING: zero testable node-split comparisons\n")
      } else {
        cat("min p:", format.pval(min(testable$p, na.rm = TRUE), digits = 4),
            " comparisons with p < 0.05:", sum(testable$p < 0.05, na.rm = TRUE), "\n")
        keep_cols <- intersect(c("comparison", "TE", "seTE", "lower", "upper",
                                 "statistic", "p", "z"), names(testable))
        print(testable[, keep_cols, drop = FALSE], row.names = FALSE, digits = 4)
      }
    }
  }

  cat("\n-- Q decomposition (design-based) --\n")
  dl <- tryCatch(netmeta::decomp.design(net), error = function(e) {
    cat("ERROR: decomp.design failed:", conditionMessage(e), "\n"); NULL
  })
  qd <- extract_qdecomp(dl)
  if (is.null(dl)) {
    cat("decomp.design unavailable\n")
  } else if (all(is.na(c(qd$Q.within, qd$Q.between)))) {
    cat("WARNING: Q decomposition fields empty; decomp names:",
        paste(names(dl), collapse = ", "), "\n")
  } else {
    cat("source=", qd$source,
        "| Q.within =", round(qd$Q.within, 2),
        "(df=", qd$df.within, ", p=",
        if (is.na(qd$p.within)) "NA" else format.pval(qd$p.within, digits = 4), ")",
        "| Q.between =", round(qd$Q.between, 2),
        "(df=", qd$df.between, ", p=",
        if (is.na(qd$p.between)) "NA" else format.pval(qd$p.between, digits = 4), ")\n")
    if (!is.null(dl$Q.decomp)) print(dl$Q.decomp, digits = 4)
  }

  if (isTRUE(save_aux) && is.null(drop_study)) {
    if (!is.null(nsdf)) {
      testable <- nsdf[!is.na(nsdf$p), , drop = FALSE]
      outfile <- if (identical(outcome_col, "responders_gi")) {
        "P2_nodesplit_gi.csv"
      } else if (identical(outcome_col, "responders_panc")) {
        "P2_nodesplit_panc.csv"
      } else {
        paste0("P2_nodesplit_", outcome_col, ".csv")
      }
      fwrite(as.data.table(testable), file.path(results_dir, outfile))
      cat("wrote", outfile, "rows:", nrow(testable), "\n")
    }
    if (!is.null(dl) && !is.null(dl$Q.decomp)) {
      qfile <- if (identical(outcome_col, "responders_gi")) {
        "P2_qdecomp_gi.csv"
      } else if (identical(outcome_col, "responders_panc")) {
        "P2_qdecomp_panc.csv"
      } else {
        paste0("P2_qdecomp_", outcome_col, ".csv")
      }
      qtab <- as.data.frame(dl$Q.decomp)
      qtab$component <- rownames(qtab)
      fwrite(as.data.table(qtab), file.path(results_dir, qfile))
      cat("wrote", qfile, "\n")
    }
  }

  invisible(list(net = net, res = res, ns = ns, nsdf = nsdf, qd = qd))
}

gi <- run_nma("responders_gi", "SERIOUS GI ADVERSE EVENTS", save_aux = TRUE)
gi_nosel <- run_nma("responders_gi", "SERIOUS GI ADVERSE EVENTS", drop_study = "NCT03574597")
pc <- run_nma("responders_panc", "PANCREATITIS", save_aux = TRUE)
pc_nosel <- run_nma("responders_panc", "PANCREATITIS", drop_study = "NCT03574597")

# save clean data + results
fwrite(DL, file.path(results_dir, "P2_clean_drug_level_data.csv"))
fwrite(gi$res, file.path(results_dir, "P2_nma_gi_corrected.csv"))
fwrite(pc$res, file.path(results_dir, "P2_nma_pancreatitis_corrected.csv"))
cat("\nDONE\n")
sink()