# -*- coding: utf-8 -*-
"""P2 NMA extractor v2: serious GI + pancreatitis with denominators, arm-label verification."""
import csv
import time
from pathlib import Path

import requests

H = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
OUT = RESULTS / "ctgov_nma_data.csv"
LOG = RESULTS / "ctgov_extract_log2.txt"

TRIALS = [
    ("NCT03548987", "semaglutide", "STEP1"), ("NCT03552757", "semaglutide", "STEP2"),
    ("NCT03611582", "semaglutide", "STEP3"), ("NCT03574597", "semaglutide", "SELECT"),
    ("NCT01720446", "semaglutide", "SUSTAIN-6"), ("NCT02692716", "semaglutide", "PIONEER-6"),
    ("NCT06131372", "semaglutide", "FLOW"),
    ("NCT01179048", "liraglutide", "LEADER"), ("NCT01272219", "liraglutide", "SCALE"),
    ("NCT03987919", "tirzepatide", "SURPASS-2"), ("NCT04184622", "tirzepatide", "SURMOUNT-1"),
    ("NCT04255433", "tirzepatide", "SURPASS-4"), ("NCT03954834", "tirzepatide", "SURMOUNT-2"),
    ("NCT01394952", "dulaglutide", "REWIND"), ("NCT01144338", "exenatide", "EXSCEL"),
    ("NCT03496298", "efpeglenatide", "AMPLITUDE-O"),
    ("NCT05931367", "retatrutide", "TRIUMPH-obesity"),
    ("NCT06297603", "retatrutide", "TRIUMPH-T2D"),
    ("NCT05869903", "orforglipron", "ATTAIN-1"),
    ("NCT06323161", "cagrisema", "REDEFINE-T2D"),
    ("NCT06066515", "survodutide", "SURVO-obesity"),
    ("NCT06164873", "mazdutide", "MAZDUTIDE-CN"),
]


def fetch(nct):
    try:
        r = requests.get(f"https://clinicaltrials.gov/api/v2/studies/{nct}", timeout=40, headers=H)
        return r.json() if r.status_code == 200 else None
    except Exception:
        return None


def extract(j):
    ps = j.get("protocolSection", {})
    nct = ps.get("identificationModule", {}).get("nctId", "?")
    title = ps.get("identificationModule", {}).get("briefTitle", "")
    rs = j.get("resultsSection", {})
    ae = rs.get("adverseEventsModule", {})
    if not ae:
        return None
    # arm labels & ids from eventGroups (module level)
    groups = {g.get("id"): g.get("title") for g in ae.get("eventGroups", [])}
    # denominators: from first serious event's stats (numAtRisk per group)
    denom = {}
    for c in ae.get("seriousEvents", []):
        for s in c.get("stats", []):
            if s.get("groupId") in groups and s.get("numAtRisk"):
                denom[s["groupId"]] = s["numAtRisk"]
    if not denom:
        return None
    gi, panc = {}, {}
    for c in ae.get("seriousEvents", []):
        st = {s["groupId"]: s.get("numAffected", 0) for s in c.get("stats", [])}
        if c.get("organSystem") == "Gastrointestinal disorders":
            for k, v in st.items():
                gi[k] = gi.get(k, 0) + v
        if "pancreatitis" in (c.get("term") or "").lower():
            for k, v in st.items():
                panc[k] = panc.get(k, 0) + v
    return {"nct": nct, "title": title, "groups": groups, "denom": denom, "gi": gi, "panc": panc}


rows = []
with open(LOG, "a", encoding="utf-8") as log:
    log.write(f"\n=== v2 {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
    for nct, drug, label in TRIALS:
        j = fetch(nct)
        if j is None:
            log.write(f"{nct}: fetch fail\n"); continue
        r = extract(j)
        if r is None:
            log.write(f"{nct}: no results/AE data\n"); continue
        title_l = r["title"].lower()
        ok = drug in title_l or any(d in title_l for d in drug.split("/"))
        log.write(f"{nct} | {label} | drug_ok={ok} | {r['title'][:70]}\n")
        for egid, t in r["groups"].items():
            # skip run-in/washout/screening groups
            if any(k in t.lower() for k in ("run-in", "washout", "screening", "lead-in")):
                continue
            n = r["denom"].get(egid)
            if n is None:
                continue
            rows.append({
                "nct": nct, "label": label, "drug_verified": ok, "arm_id": egid,
                "arm": t[:60], "arm_n": n,
                "serious_gi_affected": r["gi"].get(egid, 0),
                "pancreatitis_affected": r["panc"].get(egid, 0),
            })
        time.sleep(0.4)

with open(OUT, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["nct", "label", "drug_verified", "arm_id", "arm", "arm_n",
                                      "serious_gi_affected", "pancreatitis_affected"])
    w.writeheader()
    w.writerows(rows)
print(f"done: {len(rows)} arm-rows")
