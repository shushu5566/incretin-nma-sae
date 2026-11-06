# P2 Supplementary Code S1

Analysis scripts and derived trial-level datasets for:

> Serious gastrointestinal adverse events and pancreatitis with incretin-based therapies for obesity and type 2 diabetes: a systematic review and network meta-analysis of randomised controlled trials

**Corresponding author:** Shu-guang Li, MD (`lsg54124 [at] hotmail [dot] com`), Department of Endocrinology, The General Hospital of the Central Theater Command, PLA, Wuhan, China.

**Companion manuscript:** citation to be added on acceptance.

## Repository contents

| Path | Description |
|------|-------------|
| `scripts/` | ClinicalTrials.gov API extract/curate (Python) and netmeta NMA / figures / node-splitting (R) |
| `results/` | Derived arm-level and drug-level CSVs plus corrected analysis logs used in the manuscript |
| `LICENSE` | MIT for analysis code |
| `README.md` | This file |

Paths in scripts resolve relative to the **repository root** (parent of `scripts/`). Data are read from and written to `results/`. Internal working notes are **not** part of this release.

## Requirements (as used for the reported analyses)

- R **4.6.0** (packages may warn if built under 4.6.1) with **netmeta (>= 3.6.1)**, **meta**, **data.table** (plus **ggplot2** / **patchwork** for figures; **dplyr** / **gemtc** only for legacy `nma_*.R` scripts)
- Python **3.10+** with **requests** for optional re-extraction from ClinicalTrials.gov (network access required)

**Version note:** netmeta 3.6.x changed how `netsplit()` and `decomp.design()` expose results (`compare.random` data frame; `Q.decomp` table). The `fix_p2_*.R` scripts target that API and also write `results/P2_nodesplit_*.csv` / `results/P2_qdecomp_*.csv` for figures. Pin exact package versions in a follow-up `renv.lock` / `requirements.txt` if you need bit-for-bit reproduction beyond the manuscript tables.

## How to reproduce (high level)

1. **Preferred for manuscript numbers:** start from the derived CSVs already in `results/` and re-run the corrected NMA entry points:
   - `Rscript scripts/fix_p2_corrected_analysis.R` (writes NMA CSVs, log, node-split and Q-decomp tables)
   - then `Rscript scripts/fix_p2_figures.R` (reads those tables for Figure 3; regenerates them if missing)
   - optional: `Rscript scripts/fix_p2_nodesplit.R` / `fix_p2_nodesplit2.R`  
   Outputs overwrite or refresh files under `results/` (figures go to `results/figures/`).
2. **Optional full re-extract (network):**  
   `python3 scripts/ctgov_extract_v2.py` → `python3 scripts/ctgov_expand.py` → `python3 scripts/curate_expanded.py`  
   `ctgov_extract_v2.py` re-fetches a fixed NCT panel (byte-stable when registry posts are unchanged).  
   `ctgov_expand.py` discovers trials via API v2 with **full pagination**, server sort `LastUpdatePostDate:asc`, then a client-side **NCT ID sort**. Earlier non-determinism came from taking only the first page / first 15 title-filtered hits per drug with the API default order—not from post-freeze registry edits. Smoke the listing without the long AE download: `python3 scripts/ctgov_expand.py --discover-only`.  
   Re-runs can still diverge if ClinicalTrials.gov later adds/removes completed trials or edits AE modules; the shipped `results/*.csv` remain the manuscript freeze.
3. **Legacy scripts** `nma_combined.R`, `nma_population.R` expect intermediate tables (`nma_gi_sae_data.csv`, `nma_gi_druglevel.csv`) that are **intentionally not** shipped in this public release. Prefer the `fix_p2_*.R` pipeline above. `nma_pancreatitis.R` can use the shipped `ctgov_*.csv` files in `results/` but is superseded by `fix_p2_corrected_analysis.R` (legacy output may contain an NA treatment row; do not use for manuscript numbers).

## What is intentionally excluded

Raw database dumps, screening spreadsheets, identifier lists, protocol-registration drafts, submission correspondence, and some intermediate curation tables (`ctgov_trial_list.txt`, `nma_gi_sae_data.csv`, `nma_gi_druglevel.csv`, optional `nma_gi_expanded.csv` used only for title-based PK/DDI trial flags). Those stay on the authors' machines. `curate_expanded.py` needs `results/ctgov_trial_list.txt` (produced by `ctgov_expand.py` during a full re-extract). Without `nma_gi_expanded.csv`, the main entry still runs; title-level PK/DDI exclusion is skipped and drug-mapping exclusions still apply. The corrected analysis log states this explicitly (do not read a title-level `EXCLUDED: 0` line as "no PK studies were present").

## Study-selection arithmetic (Figure 1)

After deduplication the analysis starts from **65** unique NCT numbers. **4** trials drop out entirely because every arm fails drug mapping (insulin–incretin combinations or hepatic-impairment cohorts only): `NCT00976391`, `NCT02229227`, `NCT03940742`, `NCT05813912`. Of the remaining **61**, **26** lack a placebo arm or an incretin arm after arm-level mapping, leaving **35** trials in the NMA (65 − 4 − 26 = 35).

## Data source and licence

All input data are publicly available from [ClinicalTrials.gov](https://clinicaltrials.gov/). This repository redistributes only derived aggregate trial-level tables, not individual patient data. Analysis code is MIT (`LICENSE`).

## Citation

Cite the companion manuscript (journal citation on acceptance) and this repository URL.
