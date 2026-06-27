#!/usr/bin/env python3
"""
backfill_variant_coords.py — repair corrupt variant coordinates from cache.

The `variants` table historically stamped one gene-level dbSNP coordinate onto
every ClinVar variant of a gene, collapsing all 683 rows to a single (wrong)
locus (chr21). `src/db/populate.py` now reads per-variant GRCh38 coordinates from
the ClinVar record itself; this script applies the same corrected parse to an
ALREADY-BUILT database so it can be repaired without a full rebuild.

It reuses the exact helpers from `src/db/populate.py` so backfill and rebuild
agree. Coordinates come from cached ClinVar responses — reproducible, no network,
no hand-seeding. After running, `scripts/audit_variant_coords.py` should PASS.

Usage:
    python3 scripts/backfill_variant_coords.py [--db PATH] [--cache DIR]
"""
import argparse
import json
import os
import sys
from glob import glob

import duckdb

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS_ROOT = os.path.dirname(HERE)
sys.path.insert(0, ATLAS_ROOT)
from src.db.populate import _clinvar_grch38_loc, _clinvar_rsid  # noqa: E402

DEFAULT_DB = os.path.join(ATLAS_ROOT, "data", "processed", "als_genomic_atlas.duckdb")
DEFAULT_CACHE = os.path.join(ATLAS_ROOT, "data", "raw", "cache")


def _find_clinvar_result(obj):
    if isinstance(obj, dict):
        if "result" in obj and isinstance(obj["result"], dict) and "uids" in obj["result"]:
            return obj["result"]
        for v in obj.values():
            r = _find_clinvar_result(v)
            if r:
                return r
    return None


def collect_variant_coords(cache_dir):
    """{variant_id: (chromosome, position, hgvs, rsid)} from cached ClinVar."""
    out = {}
    files = [f for f in glob(os.path.join(cache_dir, "**", "*"), recursive=True)
             if os.path.isfile(f)]
    n_records = 0
    for f in files:
        try:
            with open(f) as fh:
                d = json.load(fh)
        except (ValueError, OSError):
            continue
        res = _find_clinvar_result(d)
        if not res:
            continue
        for uid, info in res.items():
            if uid == "uids" or not isinstance(info, dict):
                continue
            variant_id = info.get("accession") or info.get("uid")
            if not variant_id:
                continue
            chrom, pos = _clinvar_grch38_loc(info)
            hgvs = info.get("title")
            rsid = _clinvar_rsid(info)
            if chrom or pos or hgvs:
                out[variant_id] = (chrom, pos, hgvs, rsid)
                n_records += 1
    return out, n_records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    args = ap.parse_args()
    if not os.path.exists(args.db):
        sys.exit(f"ERROR: atlas DB not found: {args.db}")

    coords, n = collect_variant_coords(args.cache)
    print(f"[backfill] parsed coordinates for {len(coords)} variants from cache")
    if not coords:
        sys.exit("ERROR: no ClinVar coordinates found in cache; refusing to proceed.")

    con = duckdb.connect(args.db)
    rows = con.execute("SELECT variant_id FROM variants").fetchall()
    filled = missing = 0
    for (vid,) in rows:
        c = coords.get(vid)
        if c is None:
            missing += 1
            continue
        chrom, pos, hgvs, rsid = c
        con.execute(
            "UPDATE variants SET chromosome=?, position=?, hgvs=?, "
            "rsid=COALESCE(?, rsid) WHERE variant_id=?",
            [chrom, pos, hgvs, rsid, vid],
        )
        filled += 1
    print(f"[backfill] updated {filled}/{len(rows)} variant rows "
          f"({missing} without a cache match)")
    n_chr = con.execute("SELECT COUNT(DISTINCT chromosome) FROM variants").fetchone()[0]
    n_pos = con.execute("SELECT COUNT(DISTINCT position) FROM variants").fetchone()[0]
    print(f"[backfill] post-fix: distinct chromosome={n_chr}, distinct position={n_pos}")
    con.close()


if __name__ == "__main__":
    main()
