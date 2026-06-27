#!/usr/bin/env python3
"""
audit_variant_coords.py — integrity gate for the `variants` coordinate columns.

WHY
---
The `variants` table's coordinate fields (`chromosome`, `position`, `hgvs`) were
found to be corrupt PANEL-WIDE: every row carries the identical placeholder
`chromosome=21`, one `position`, and `hgvs='NC_000021.9:g.31659783C>T'`,
regardless of the gene's real locus. Any downstream use of variant coordinates
(coordinate-based annotation, regulatory overlap, plotting) is therefore
meaningless. The per-gene scalar annotations (gnomad_pli/loeuf) are unaffected.

This script is a CI-style GATE: it cross-checks variant coordinates against the
gene's real chromosome (from the `genes` table) and against the expectation that
distinct variants have distinct positions, then exits non-zero if corruption is
detected so a broken rebuild fails loudly instead of shipping bad coordinates.

It does NOT mutate the DB. The correct fix lives upstream in the ClinVar variant
ingestion (the coordinate parse is writing a constant); until that is fixed,
variant coordinates should be treated as UNTRUSTWORTHY and re-derived from the
live ClinVar API per variant.

Usage:
    python3 scripts/audit_variant_coords.py [--db PATH]
Exit code 0 = coordinates look sane; 1 = corruption detected.
"""
import argparse
import os
import sys

import duckdb

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS_ROOT = os.path.dirname(HERE)
DEFAULT_DB = os.path.join(ATLAS_ROOT, "data", "processed", "als_genomic_atlas.duckdb")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=DEFAULT_DB)
    args = ap.parse_args()
    if not os.path.exists(args.db):
        sys.exit(f"ERROR: atlas DB not found: {args.db}")

    con = duckdb.connect(args.db, read_only=True)

    total = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    n_chr = con.execute("SELECT COUNT(DISTINCT chromosome) FROM variants").fetchone()[0]
    n_pos = con.execute("SELECT COUNT(DISTINCT position) FROM variants").fetchone()[0]
    n_hgvs = con.execute("SELECT COUNT(DISTINCT hgvs) FROM variants").fetchone()[0]
    n_genes_with_vars = con.execute(
        "SELECT COUNT(DISTINCT gene_symbol) FROM variants").fetchone()[0]

    print(f"[audit] variants={total} across {n_genes_with_vars} genes | "
          f"distinct chromosome={n_chr} position={n_pos} hgvs={n_hgvs}")

    problems = []

    # Signal 1: collapsed coordinates (a healthy table has many distinct
    # positions; one shared position/hgvs across all genes is the failure mode).
    if total > 0 and n_pos <= 1:
        problems.append(f"all {total} variants share a single `position` "
                        f"(distinct positions = {n_pos}).")
    if total > 0 and n_hgvs <= 1:
        problems.append(f"all {total} variants share a single `hgvs` value.")
    if n_genes_with_vars > 1 and n_chr <= 1:
        problems.append(f"variants span {n_genes_with_vars} genes but only "
                        f"{n_chr} distinct chromosome(s) — genes are on many chromosomes.")

    # Signal 2: per-gene chromosome mismatch vs the authoritative `genes` table.
    mism = con.execute("""
        SELECT g.gene_symbol, g.chromosome AS gene_chr,
               COUNT(*) AS nvar,
               COUNT(DISTINCT v.chromosome) AS var_chr_vals,
               MIN(v.chromosome) AS sample_var_chr
        FROM genes g
        JOIN variants v ON v.gene_symbol = g.gene_symbol
        GROUP BY 1, 2
        HAVING MIN(v.chromosome) IS DISTINCT FROM CAST(g.chromosome AS VARCHAR)
        ORDER BY 1
    """).fetchall()
    if mism:
        problems.append(f"{len(mism)} genes have variants whose chromosome "
                        "does not match the gene's chromosome (sample below).")

    if mism:
        print("\n[audit] gene-chromosome mismatches (first 10):")
        print(f"  {'gene':<10}{'gene_chr':<10}{'var_chr':<10}{'nvar':<6}")
        for gene, gchr, nvar, _vals, vchr in mism[:10]:
            print(f"  {gene:<10}{str(gchr):<10}{str(vchr):<10}{nvar:<6}")

    con.close()

    if problems:
        print("\n[audit] FAIL — variant coordinates are CORRUPT:")
        for p in problems:
            print(f"  - {p}")
        print("\n  Fix upstream in the ClinVar variant ingestion (coordinate "
              "parse is writing a constant). Until then, treat variant "
              "chromosome/position/hgvs as UNTRUSTWORTHY and re-derive from the "
              "live ClinVar API. (Per-gene gnomad_pli/loeuf are unaffected.)")
        sys.exit(1)

    print("[audit] PASS — variant coordinates look sane.")
    sys.exit(0)


if __name__ == "__main__":
    main()
