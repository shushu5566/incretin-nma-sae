# -*- coding: utf-8 -*-
"""P2: curate expanded CT.gov data -> NMA-ready dataset (serious GI SAE)."""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DIR = RESULTS
titles = {}
with open(RESULTS / "ctgov_trial_list.txt", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            titles[p[0]] = p[1]

with open(RESULTS / "ctgov_expanded_data.csv", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

# ---- trial-level exclusion ----
EXCLUDE_TITLE = ["drug-drug interaction", "pharmacokinetic", "device", "healthy participant",
                 "healthy volunteer", "impaired liver", "impaired renal", "qt interval",
                 "time in range", "switching", "icodec", "rdx-002", "type 1 diabetes",
                 "dose ranging", "supratherapeutic", "master protocol"]
# RISE = pediatric study (NCT01779362 - RISE Adult Medication Study - actually adult; keep? RISE is T2D risk adults - keep)
keep_nct = set()
for r in rows:
    nct = r["nct"]
    t = titles.get(nct, "").lower()
    if any(k in t for k in EXCLUDE_TITLE):
        continue
    keep_nct.add(nct)

# ---- arm -> node mapping (drug x dose) ----
def map_node(arm):
    a = arm.lower()
    # placebo variants
    if "placebo" in a:
        return "Placebo"
    # semaglutide
    m = re.search(r"semaglutide\s*([0-9.]+)\s*mg", a)
    if m:
        return f"Sema{m.group(1)}"
    if "semaglutide" in a and "oral" in a:
        return "OralSema"
    if "semaglutide" in a:
        return "Sema"
    # liraglutide
    m = re.search(r"liraglutide\s*([0-9.]+)\s*mg", a)
    if m:
        return f"Lira{m.group(1)}"
    if "liraglutide" in a:
        return "Lira"
    # tirzepatide (LY3298176)
    m = re.search(r"tirzepatide[^0-9]*([0-9.]+)\s*mg", a)
    if m:
        return f"Tirz{m.group(1)}"
    if "tirzepatide" in a or "ly3298176" in a:
        return "Tirz"
    # dulaglutide (LY2189265)
    m = re.search(r"dulaglutide[^0-9]*([0-9.]+)\s*mg", a)
    if m:
        return f"Dula{m.group(1)}"
    if "dulaglutide" in a or "ly2189265" in a:
        return "Dula"
    # exenatide
    m = re.search(r"exenatide[^0-9]*([0-9.]+)\s*(?:mg|mcg|ug)", a)
    if m:
        return f"Exe{m.group(1)}"
    if "exenatide" in a:
        return "Exe"
    # lixisenatide
    if "lixisenatide" in a:
        return "Lixi"
    # albiglutide
    if "albiglutide" in a:
        return "Albig"
    # efpeglenatide
    if "efpeglenatide" in a:
        return "Efpe"
    # orforglipron (LY3502970)
    m = re.search(r"(?:ly3502970|orforglipron)[^0-9]*([0-9.]+)\s*mg", a)
    if m:
        return f"Orfor{m.group(1)}"
    if "ly3502970" in a or "orforglipron" in a:
        return "Orfor"
    # retatrutide (LY3437943)
    m = re.search(r"(?:ly3437943|retatrutide)[^0-9]*([0-9.]+)\s*mg", a)
    if m:
        return f"Reta{m.group(1)}"
    if "ly3437943" in a or "retatrutide" in a:
        return "Reta"
    # survodutide (BI 456906)
    if "456906" in a or "survodutide" in a:
        return "Survo"
    # sitagliptin (non-incretin comparator - exclude)
    return None

out = []
skip = 0
for r in rows:
    if r["nct"] not in keep_nct:
        continue
    node = map_node(r["arm"])
    if node is None:
        skip += 1
        continue
    out.append({"nct": r["nct"], "trial_title": titles.get(r["nct"], ""), "arm": r["arm"],
                "node": node, "n": int(r["arm_n"]), "gi": int(r["serious_gi_affected"])})

# pool identical (trial, node)
from collections import defaultdict
agg = defaultdict(lambda: {"n": 0, "gi": 0})
for o in out:
    k = (o["nct"], o["node"])
    agg[k]["n"] += o["n"]
    agg[k]["gi"] += o["gi"]

# keep trials with >=2 arms after pooling
by_trial = defaultdict(list)
for (nct, node), v in agg.items():
    by_trial[nct].append((node, v["n"], v["gi"]))
final = []
for nct, arms in by_trial.items():
    if len(arms) < 2:
        continue
    for node, n, gi in arms:
        final.append({"study": nct, "treatment": node, "responders": gi, "sampleSize": n,
                      "title": titles.get(nct, "")})

with open(RESULTS / "nma_gi_expanded.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["study", "treatment", "responders", "sampleSize", "title"])
    w.writeheader()
    w.writerows(final)
print("mapped rows:", len(out), "| unmapped arms:", skip)
print("final NMA arms:", len(final), "| trials:", len(set(x["study"] for x in final)),
      "| nodes:", len(set(x["treatment"] for x in final)))
from collections import Counter
print("node counts:", dict(Counter(x["treatment"] for x in final)))
