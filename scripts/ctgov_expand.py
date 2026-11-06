# -*- coding: utf-8 -*-
"""P2: discover phase-3 incretin trials on CT.gov and batch-extract AE data.

Discovery paginates ClinicalTrials.gov API v2 until exhausted, requests a stable
server-side sort (LastUpdatePostDate:asc), then re-sorts hits by NCT ID client-side.
Use --discover-only to smoke-test listing without the long AE download.
"""
import argparse
import csv
import time
from pathlib import Path

import requests

H = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
DIR = RESULTS
LOG = RESULTS / "ctgov_expand_log.txt"

DRUGS = ["semaglutide", "liraglutide", "tirzepatide", "dulaglutide", "exenatide",
         "lixisenatide", "albiglutide", "efpeglenatide", "orforglipron",
         "retatrutide", "survodutide", "mazdutide", "cagrilintide"]

# Stable server-side sort so page tokens do not reshuffle mid-pagination.
# After all pages are collected, hits are also sorted by NCT ID client-side.
API_SORT = "LastUpdatePostDate:asc"
PAGE_SIZE = 100
TITLE_KEYS = ("trial", "study", "effect", "compare", "safety", "efficacy",
              "cardiovascular", "obesity", "diabetes")


def discover(drug, maxn=None):
    """Return [(nctId, briefTitle), ...] for one intervention term.

    Paginates with nextPageToken (no first-page truncation). Optional maxn is
    applied only AFTER full collection + NCT-ID sort (for dry-run smoke tests).
    """
    ids = []
    page_token = None
    while True:
        params = {
            "query.intr": drug,
            "filter.overallStatus": "COMPLETED",
            "pageSize": PAGE_SIZE,
            "sort": API_SORT,
            "fields": "NCTId,BriefTitle,StudyType",
        }
        if page_token:
            params["pageToken"] = page_token
        try:
            r = requests.get("https://clinicaltrials.gov/api/v2/studies",
                             params=params, timeout=60, headers=H)
            if r.status_code != 200:
                with open(LOG, "a", encoding="utf-8") as f:
                    f.write(f"discover {drug} HTTP {r.status_code}\\n")
                break
            data = r.json()
        except Exception as e:
            with open(LOG, "a", encoding="utf-8") as f:
                f.write(f"discover {drug} FAIL {str(e)[:80]}\\n")
            break
        for s in data.get("studies", []):
            ps = s.get("protocolSection", {})
            ident = ps.get("identificationModule", {})
            title = (ident.get("briefTitle", "") or "").lower()
            stype = ps.get("designModule", {}).get("studyType", "")
            if stype == "INTERVENTIONAL" and any(k in title for k in TITLE_KEYS):
                ids.append((ident["nctId"], ident.get("briefTitle", "")[:90]))
        page_token = data.get("nextPageToken")
        if not page_token:
            break
        time.sleep(0.25)
    # Client-side stable order independent of API default ranking
    ids = sorted(set(ids), key=lambda x: x[0])
    if maxn is not None:
        ids = ids[:maxn]
    return ids


def fetch(nct):
    try:
        r = requests.get(f"https://clinicaltrials.gov/api/v2/studies/{nct}",
                         timeout=40, headers=H)
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
    groups = {g.get("id"): g.get("title") for g in ae.get("eventGroups", [])}
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--discover-only", action="store_true",
                    help="List unique NCT IDs then exit (no AE extraction)")
    ap.add_argument("--maxn-per-drug", type=int, default=None,
                    help="Optional per-drug cap AFTER full pagination+sort (smoke tests)")
    args = ap.parse_args()

    with open(LOG, "a", encoding="utf-8") as log:
        log.write(f"\\n=== discovery {time.strftime('%Y-%m-%d %H:%M:%S')} ===\\n")
        log.write(f"sort={API_SORT} pageSize={PAGE_SIZE} maxn_per_drug={args.maxn_per_drug}\\n")
        all_trials = []
        for drug in DRUGS:
            hits = discover(drug, maxn=args.maxn_per_drug)
            all_trials.extend(hits)
            log.write(f"{drug}: {len(hits)} trials\\n")
            time.sleep(0.4)
        seen = set()
        uniq = []
        for nct, t in all_trials:
            if nct not in seen:
                seen.add(nct)
                uniq.append((nct, t))
        uniq.sort(key=lambda x: x[0])
        log.write(f"total unique: {len(uniq)}\\n")
        with open(RESULTS / "ctgov_trial_list.txt", "w", encoding="utf-8") as f:
            for nct, t in uniq:
                f.write(f"{nct}\\t{t}\\n")
        print(f"discover: {len(uniq)} unique trials (list written)")
        if args.discover_only:
            log.write("discover-only: skip AE extraction\\n")
            print("discover-only: done")
            return

        rows = []
        log.write("=== extraction ===\\n")
        for nct, t in uniq:
            j = fetch(nct)
            if j is None:
                log.write(f"{nct} fetch fail\\n")
                continue
            r = extract(j)
            if r is None:
                log.write(f"{nct} no AE results\\n")
                continue
            log.write(f"{nct} | {r['title'][:70]}\\n")
            for egid, atitle in r["groups"].items():
                if any(k in atitle.lower() for k in ("run-in", "washout", "screening", "lead-in")):
                    continue
                n = r["denom"].get(egid)
                if n is None:
                    continue
                rows.append({"nct": nct, "arm_id": egid, "arm": atitle[:60], "arm_n": n,
                             "serious_gi_affected": r["gi"].get(egid, 0),
                             "pancreatitis_affected": r["panc"].get(egid, 0)})
            time.sleep(0.3)
        log.write(f"arm rows: {len(rows)}\\n")
        with open(RESULTS / "ctgov_expanded_data.csv", "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["nct", "arm_id", "arm", "arm_n",
                                              "serious_gi_affected", "pancreatitis_affected"])
            w.writeheader()
            w.writerows(rows)
    print(f"done: {len(rows)} arm rows from {len(uniq)} trials")


if __name__ == "__main__":
    main()
