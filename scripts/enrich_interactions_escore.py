#!/usr/bin/env python3
"""
enrich_interactions_escore.py — backfill STRING evidence channels into the atlas.

WHY
---
`populate_interactions` historically stored only the STRING *combined* score in
`interactions.confidence_score`. The combined score conflates evidence channels:
an edge can score ~0.9 purely from TEXT-MINING co-mentions (escore~0, tscore~1)
with NO physical/functional evidence. Cross-dossier verification found that of
the panel-gene "convergences" implied by the combined score, only ~5 are genuine
(experimental/curated) and ~8 are text-mining artifacts.

The per-edge sub-scores (escore=experimental, dscore=database/curated,
tscore=text-mining) are ALREADY present in the cached STRING responses — they
were simply discarded at populate time. This script:

  1. Reads the STRING evidence channels back out of `data/raw/cache` (no network,
     fully reproducible — same source `populate` consumes).
  2. Adds escore/dscore/tscore columns to the existing `interactions` table and
     backfills them (idempotent; safe to re-run). New full rebuilds get the
     columns natively via the updated schema.py/populate.py.
  3. Writes a panel<->panel "genuine vs text-mining convergence" report so the
     misleading edges are explicit.

Reproducibility: every value comes from a cached STRING API response; nothing is
hand-seeded. Run against the same cache twice -> identical output.

Usage:
    python3 scripts/enrich_interactions_escore.py [--db PATH] [--cache DIR] [--report PATH]
"""
import argparse
import json
import os
import sys
from glob import glob

import duckdb

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS_ROOT = os.path.dirname(HERE)
DEFAULT_DB = os.path.join(ATLAS_ROOT, "data", "processed", "als_genomic_atlas.duckdb")
DEFAULT_CACHE = os.path.join(ATLAS_ROOT, "data", "raw", "cache")
DEFAULT_REPORT = os.path.join(ATLAS_ROOT, "outputs", "interaction_provenance_report.md")

# Classification thresholds (STRING channel conventions).
ESCORE_GENUINE = 0.40   # experimental support at/above this = physically supported
DSCORE_CURATED = 0.70   # curated database support at/above this = trustworthy
TSCORE_TEXTMINING = 0.50  # text-mining at/above this with low escore = co-mention


def _find_string_edges(cache_dir):
    """Walk the cache and collect STRING interaction edges with sub-scores.

    Returns {(gene_a, gene_b sorted): {score, escore, dscore, tscore}}.
    """
    edges = {}
    files = [f for f in glob(os.path.join(cache_dir, "**", "*"), recursive=True)
             if os.path.isfile(f)]

    def find_string_list(obj):
        if isinstance(obj, list) and obj and isinstance(obj[0], dict) \
                and "preferredName_A" in obj[0] and "escore" in obj[0]:
            return obj
        if isinstance(obj, dict):
            for v in obj.values():
                r = find_string_list(v)
                if r:
                    return r
        return None

    n_responses = 0
    for f in files:
        try:
            with open(f) as fh:
                d = json.load(fh)
        except (ValueError, OSError):
            continue
        sl = find_string_list(d)
        if not sl:
            continue
        n_responses += 1
        for item in sl:
            a, b = item.get("preferredName_A"), item.get("preferredName_B")
            if not (a and b):
                continue
            key = tuple(sorted([a, b]))
            # STRING is symmetric; if seen twice (A's and B's query) values match.
            edges[key] = {
                "score": item.get("score"),
                "escore": item.get("escore"),
                "dscore": item.get("dscore"),
                "tscore": item.get("tscore"),
            }
    return edges, n_responses


def _classify(escore, dscore, tscore):
    e = escore or 0.0
    d = dscore or 0.0
    t = tscore or 0.0
    if e >= ESCORE_GENUINE or d >= DSCORE_CURATED:
        return "GENUINE"
    if e < 0.10 and t >= TSCORE_TEXTMINING:
        return "TEXT-MINING"
    return "WEAK"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    args = ap.parse_args()

    if not os.path.exists(args.db):
        sys.exit(f"ERROR: atlas DB not found: {args.db}")
    if not os.path.isdir(args.cache):
        sys.exit(f"ERROR: cache dir not found: {args.cache}")

    edges, n_responses = _find_string_edges(args.cache)
    print(f"[enrich] read {len(edges)} STRING edges from {n_responses} cached responses")
    if not edges:
        sys.exit("ERROR: no STRING edges found in cache; refusing to proceed "
                 "(would otherwise null out the columns).")

    con = duckdb.connect(args.db)  # writable
    # Idempotent column adds.
    for col in ("escore", "dscore", "tscore"):
        con.execute(f"ALTER TABLE interactions ADD COLUMN IF NOT EXISTS {col} DOUBLE")

    # Backfill from cache by sorted pair.
    rows = con.execute("SELECT gene_a, gene_b FROM interactions").fetchall()
    filled = 0
    missing = 0
    for a, b in rows:
        key = tuple(sorted([a, b]))
        sub = edges.get(key)
        if sub is None:
            missing += 1
            continue
        con.execute(
            "UPDATE interactions SET escore=?, dscore=?, tscore=? WHERE gene_a=? AND gene_b=?",
            [sub["escore"], sub["dscore"], sub["tscore"], a, b],
        )
        filled += 1
    print(f"[enrich] backfilled {filled}/{len(rows)} interaction rows "
          f"({missing} without a cache match)")

    # Panel<->panel convergence report (both endpoints are seed genes).
    panel = {r[0] for r in con.execute("SELECT gene_symbol FROM genes").fetchall()}
    pp = con.execute("""
        SELECT gene_a, gene_b, confidence_score, escore, dscore, tscore
        FROM interactions
        WHERE gene_a IN (SELECT gene_symbol FROM genes)
          AND gene_b IN (SELECT gene_symbol FROM genes)
        ORDER BY escore DESC NULLS LAST, confidence_score DESC
    """).fetchall()

    classified = []
    for a, b, sc, e, d, t in pp:
        classified.append((a, b, sc, e, d, t, _classify(e, d, t)))

    counts = {"GENUINE": 0, "WEAK": 0, "TEXT-MINING": 0}
    for *_, cls in classified:
        counts[cls] += 1

    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, "w") as fh:
        fh.write("# Atlas interaction-provenance report (panel <-> panel edges)\n\n")
        fh.write("Generated by `scripts/enrich_interactions_escore.py` from cached "
                 "STRING responses. Classification uses STRING evidence channels, "
                 "NOT the combined score:\n\n")
        fh.write(f"- **GENUINE**: escore >= {ESCORE_GENUINE} (experimental) OR "
                 f"dscore >= {DSCORE_CURATED} (curated)\n")
        fh.write(f"- **TEXT-MINING**: escore < 0.10 AND tscore >= {TSCORE_TEXTMINING} "
                 "(co-mention only; NOT physical evidence)\n")
        fh.write("- **WEAK**: everything else (sub-threshold support)\n\n")
        fh.write(f"**Summary:** {counts['GENUINE']} GENUINE, {counts['WEAK']} WEAK, "
                 f"{counts['TEXT-MINING']} TEXT-MINING "
                 f"(of {len(classified)} panel<->panel edges).\n\n")
        fh.write("| gene_a | gene_b | combined | escore | dscore | tscore | class |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for a, b, sc, e, d, t, cls in classified:
            fh.write(f"| {a} | {b} | {fmt(sc)} | {fmt(e)} | {fmt(d)} | {fmt(t)} | "
                     f"{'**'+cls+'**' if cls!='WEAK' else cls} |\n")
        fh.write("\n> A `confidence_score` near 1.0 with `escore`≈0 and `tscore`≈1 is a "
                 "literature co-mention, not a physical/functional interaction. Filter on "
                 "`escore`/`dscore` before asserting any panel convergence.\n")

    print(f"[enrich] wrote convergence report -> {args.report}")
    print(f"[enrich] panel<->panel: GENUINE={counts['GENUINE']} "
          f"WEAK={counts['WEAK']} TEXT-MINING={counts['TEXT-MINING']}")
    con.close()


def fmt(x):
    return "—" if x is None else f"{x:.3f}"


if __name__ == "__main__":
    main()
