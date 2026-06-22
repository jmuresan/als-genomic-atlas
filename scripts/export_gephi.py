#!/usr/bin/env python3
"""
Export DuckDB networks to Gephi-friendly CSV pairs (nodes + edges).

One-shot import (recommended):
  python scripts/export_gephi.py --network unified

  -> outputs/gephi/atlas_unified_edges.csv   (column `network` tags each layer)
  -> outputs/gephi/atlas_unified_nodes.csv   (shared node ids; als_seed=1 on panel genes)

Layers in unified graph: string, gene-pathway, foldseek, gene-bodyimpact, gene-disease

Per-layer CSVs: --network all --split   (or individual network names)
"""

from __future__ import annotations

import argparse
import os
import sys

import duckdb

DEFAULT_DB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "processed",
    "als_genomic_atlas.duckdb",
)
DEFAULT_OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "outputs",
    "gephi",
)

PANEL_DISEASE_LABEL = "Amyotrophic lateral sclerosis"
PANEL_DISEASE_ID = f"DISEASE:{PANEL_DISEASE_LABEL}"


def _export_csv(conn: duckdb.DuckDBPyConnection, sql: str, path: str) -> int:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    conn.execute(
        f"COPY ({sql}) TO ? (HEADER, DELIMITER ',')",
        [path],
    )
    row = conn.execute(f"SELECT COUNT(*) FROM ({sql}) t").fetchone()
    return int(row[0]) if row else 0


def _sql_normalize_by_network(raw_sql: str, trailing_select: str = "") -> str:
    """Min–max scale Weight to [0, 1] separately for each `network` value.
    If all raw weights in a layer are equal, Weight is set to 0.33."""
    trail = f", {trailing_select.lstrip(', ')}" if trailing_select.strip() else ""
    return f"""
    WITH raw_edges AS (
        {raw_sql}
    ),
    scaled AS (
        SELECT
            *,
            min(Weight) OVER (PARTITION BY network) AS w_min,
            max(Weight) OVER (PARTITION BY network) AS w_max
        FROM raw_edges
    )
    SELECT
        Source,
        Target,
        Weight AS WeightRaw,
        CASE
            WHEN w_max = w_min THEN 0.33
            ELSE (Weight - w_min) / (w_max - w_min)
        END AS Weight,
        network,
        EdgeType{trail}
    FROM scaled
    """


def _foldseek_where(min_probability: float | None) -> str:
    if min_probability is None:
        return ""
    return f"WHERE probability >= {float(min_probability)}"


def export_gene_pathway(conn: duckdb.DuckDBPyConnection, out_dir: str) -> None:
    edges_sql = _sql_normalize_by_network(
        """
        SELECT
            gene_symbol AS Source,
            id AS Target,
            1.0 AS Weight,
            'gene-pathway' AS network,
            type AS EdgeType,
            name AS TargetName
        FROM pathways_and_domains
        """,
        "TargetName",
    )
    nodes_sql = """
        WITH genes AS (
            SELECT DISTINCT gene_symbol AS Id FROM pathways_and_domains
        ),
        terms AS (
            SELECT DISTINCT
                id AS Id,
                name AS Label,
                type AS node_type
            FROM pathways_and_domains
        )
        SELECT Id, Id AS Label, 'gene' AS node_type, NULL AS term_type, NULL AS term_name
        FROM genes
        UNION ALL
        SELECT Id, COALESCE(Label, Id) AS Label, 'annotation' AS node_type, node_type AS term_type, Label AS term_name
        FROM terms
    """
    e_path = os.path.join(out_dir, "gene_pathway_edges.csv")
    n_path = os.path.join(out_dir, "gene_pathway_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"gene-pathway: {ne} edges -> {e_path}")
    print(f"gene-pathway: {nn} nodes -> {n_path}")


def export_string(conn: duckdb.DuckDBPyConnection, out_dir: str) -> None:
    edges_sql = _sql_normalize_by_network(
        """
        SELECT
            least(gene_a, gene_b) AS Source,
            greatest(gene_a, gene_b) AS Target,
            max(confidence_score) AS Weight,
            'string' AS network,
            'ppi' AS EdgeType
        FROM interactions
        GROUP BY 1, 2
        """
    )
    nodes_sql = """
        WITH edge_nodes AS (
            SELECT gene_a AS Id FROM interactions
            UNION
            SELECT gene_b FROM interactions
        )
        SELECT
            n.Id,
            n.Id AS Label,
            CASE WHEN g.gene_symbol IS NOT NULL THEN 'als_panel_gene' ELSE 'string_partner_gene' END AS node_type,
            CASE WHEN g.gene_symbol IS NOT NULL THEN 1 ELSE 0 END AS als_seed
        FROM edge_nodes n
        LEFT JOIN genes g ON g.gene_symbol = n.Id
    """
    e_path = os.path.join(out_dir, "string_edges.csv")
    n_path = os.path.join(out_dir, "string_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"string: {ne} edges -> {e_path}")
    print(f"string: {nn} nodes -> {n_path}")


def export_foldseek(
    conn: duckdb.DuckDBPyConnection,
    out_dir: str,
    min_probability: float | None,
) -> None:
    where = _foldseek_where(min_probability)
    edges_sql = _sql_normalize_by_network(
        f"""
        SELECT
            query_gene_symbol AS Source,
            target_id AS Target,
            probability AS Weight,
            'foldseek' AS network,
            db AS EdgeType,
            query_coverage AS QueryCoverage,
            evalue AS EValue,
            seq_identity AS SeqIdentity,
            alignment_length AS AlnLength
        FROM foldseek_matches
        {where}
        """,
        "QueryCoverage, EValue, SeqIdentity, AlnLength",
    )
    nodes_sql = f"""
        WITH q AS (
            SELECT DISTINCT query_gene_symbol AS Id FROM foldseek_matches {where}
        ),
        t AS (
            SELECT DISTINCT target_id AS Id FROM foldseek_matches {where}
        )
        SELECT Id, Id AS Label, 'query_gene' AS node_type FROM q
        UNION ALL
        SELECT
            Id,
            COALESCE(
                NULLIF(regexp_extract(Id, '(AF-[A-Z0-9]+)', 1), ''),
                LEFT(Id, 80)
            ) AS Label,
            'foldseek_target' AS node_type
        FROM t
    """
    e_path = os.path.join(out_dir, "foldseek_edges.csv")
    n_path = os.path.join(out_dir, "foldseek_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"foldseek: {ne} edges -> {e_path}")
    print(f"foldseek: {nn} nodes -> {n_path}")


def export_gene_bodyimpact(conn: duckdb.DuckDBPyConnection, out_dir: str) -> None:
    edges_sql = _sql_normalize_by_network(
        """
        SELECT
            gene_symbol AS Source,
            tissue AS Target,
            tpm AS Weight,
            'gene-bodyimpact' AS network,
            'expression' AS EdgeType
        FROM expression
        WHERE tpm IS NOT NULL
        """
    )
    nodes_sql = """
        WITH g AS (
            SELECT DISTINCT gene_symbol AS Id FROM expression
        ),
        t AS (
            SELECT DISTINCT tissue AS Id FROM expression
        )
        SELECT Id, Id AS Label, 'gene' AS node_type FROM g
        UNION ALL
        SELECT Id, Id AS Label, 'tissue' AS node_type FROM t
    """
    e_path = os.path.join(out_dir, "gene_bodyimpact_edges.csv")
    n_path = os.path.join(out_dir, "gene_bodyimpact_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"gene-bodyimpact: {ne} edges -> {e_path}")
    print(f"gene-bodyimpact: {nn} nodes -> {n_path}")


def export_gene_disease(conn: duckdb.DuckDBPyConnection, out_dir: str) -> None:
    edges_sql = _sql_normalize_by_network(
        f"""
        SELECT Source, Target, Weight, network, EdgeType FROM (
            SELECT
                gene_symbol AS Source,
                'DISEASE:' || disease_name AS Target,
                COUNT(*)::DOUBLE AS Weight,
                'gene-disease' AS network,
                'clinvar' AS EdgeType
            FROM variants
            WHERE disease_name IS NOT NULL AND trim(disease_name) <> ''
            GROUP BY gene_symbol, disease_name

            UNION ALL

            SELECT
                g.gene_symbol AS Source,
                '{PANEL_DISEASE_ID}' AS Target,
                1.0 AS Weight,
                'gene-disease' AS network,
                'seed_panel' AS EdgeType
            FROM genes g
            WHERE NOT EXISTS (
                SELECT 1 FROM variants v WHERE v.gene_symbol = g.gene_symbol
            )
        ) u
        """
    )
    nodes_sql = f"""
        WITH gene_nodes AS (
            SELECT
                g.gene_symbol AS Id,
                g.gene_symbol AS Label,
                'gene' AS node_type,
                1 AS als_seed,
                g.ensembl_id,
                g.uniprot_id,
                g.protein_description,
                COALESCE(vc.n, 0)::INTEGER AS clinvar_row_count
            FROM genes g
            LEFT JOIN (
                SELECT gene_symbol, COUNT(*) AS n FROM variants GROUP BY gene_symbol
            ) vc ON vc.gene_symbol = g.gene_symbol
        ),
        disease_nodes AS (
            SELECT DISTINCT
                'DISEASE:' || disease_name AS Id,
                disease_name AS Label,
                'disease' AS node_type,
                0 AS als_seed,
                NULL AS ensembl_id,
                NULL AS uniprot_id,
                NULL AS protein_description,
                NULL::INTEGER AS clinvar_row_count
            FROM variants
            WHERE disease_name IS NOT NULL AND trim(disease_name) <> ''
        ),
        panel_disease AS (
            SELECT
                '{PANEL_DISEASE_ID}' AS Id,
                '{PANEL_DISEASE_LABEL}' AS Label,
                'disease' AS node_type,
                0 AS als_seed,
                NULL AS ensembl_id,
                NULL AS uniprot_id,
                'ALS seed gene panel (reference disease)' AS protein_description,
                NULL::INTEGER AS clinvar_row_count
        )
        SELECT * FROM gene_nodes
        UNION ALL
        SELECT * FROM disease_nodes
        UNION ALL
        SELECT * FROM panel_disease
        WHERE NOT EXISTS (SELECT 1 FROM disease_nodes d WHERE d.Id = '{PANEL_DISEASE_ID}')
    """
    e_path = os.path.join(out_dir, "gene_disease_edges.csv")
    n_path = os.path.join(out_dir, "gene_disease_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"gene-disease: {ne} edges -> {e_path}")
    print(f"gene-disease: {nn} nodes -> {n_path}")


def export_unified(
    conn: duckdb.DuckDBPyConnection,
    out_dir: str,
    min_probability: float | None,
    include_foldseek: bool,
) -> None:
    """Single graph: shared gene nodes link all layers via the `network` edge column."""
    fs_where = _foldseek_where(min_probability)
    foldseek_edges = ""
    foldseek_nodes = ""
    if include_foldseek:
        foldseek_edges = f"""
        UNION ALL
        SELECT
            query_gene_symbol AS Source,
            target_id AS Target,
            probability AS Weight,
            'foldseek' AS network,
            db AS EdgeType
        FROM foldseek_matches
        {fs_where}
        """
        foldseek_nodes = f"""
        UNION ALL
        SELECT
            target_id AS Id,
            COALESCE(
                NULLIF(regexp_extract(target_id, '(AF-[A-Z0-9]+)', 1), ''),
                LEFT(target_id, 80)
            ) AS Label,
            'foldseek_target' AS node_type,
            0 AS als_seed,
            NULL AS ensembl_id,
            NULL AS uniprot_id,
            NULL AS protein_description,
            NULL::INTEGER AS clinvar_row_count
        FROM (SELECT DISTINCT target_id FROM foldseek_matches {fs_where}) t
        """

    edges_sql = _sql_normalize_by_network(
        f"""
        SELECT
            least(gene_a, gene_b) AS Source,
            greatest(gene_a, gene_b) AS Target,
            max(confidence_score) AS Weight,
            'string' AS network,
            'ppi' AS EdgeType
        FROM interactions
        GROUP BY 1, 2

        UNION ALL

        SELECT
            gene_symbol AS Source,
            id AS Target,
            1.0 AS Weight,
            'gene-pathway' AS network,
            type AS EdgeType
        FROM pathways_and_domains

        UNION ALL

        SELECT
            gene_symbol AS Source,
            tissue AS Target,
            tpm AS Weight,
            'gene-bodyimpact' AS network,
            'expression' AS EdgeType
        FROM expression
        WHERE tpm IS NOT NULL

        UNION ALL

        SELECT
            gene_symbol AS Source,
            'DISEASE:' || disease_name AS Target,
            COUNT(*)::DOUBLE AS Weight,
            'gene-disease' AS network,
            'clinvar' AS EdgeType
        FROM variants
        WHERE disease_name IS NOT NULL AND trim(disease_name) <> ''
        GROUP BY gene_symbol, disease_name

        UNION ALL

        SELECT
            g.gene_symbol AS Source,
            '{PANEL_DISEASE_ID}' AS Target,
            1.0 AS Weight,
            'gene-disease' AS network,
            'seed_panel' AS EdgeType
        FROM genes g
        WHERE NOT EXISTS (SELECT 1 FROM variants v WHERE v.gene_symbol = g.gene_symbol)

        {foldseek_edges}
        """
    )

    nodes_sql = f"""
        WITH combined AS (
            SELECT
                g.gene_symbol AS Id,
                g.gene_symbol AS Label,
                'gene' AS node_type,
                1 AS als_seed,
                g.ensembl_id,
                g.uniprot_id,
                g.protein_description,
                COALESCE(vc.n, 0)::INTEGER AS clinvar_row_count
            FROM genes g
            LEFT JOIN (
                SELECT gene_symbol, COUNT(*) AS n FROM variants GROUP BY gene_symbol
            ) vc ON vc.gene_symbol = g.gene_symbol

            UNION ALL

            SELECT
                n.Id,
                n.Id AS Label,
                'gene' AS node_type,
                0 AS als_seed,
                NULL, NULL, NULL, NULL::INTEGER
            FROM (
                SELECT gene_a AS Id FROM interactions
                UNION
                SELECT gene_b FROM interactions
            ) n
            WHERE n.Id NOT IN (SELECT gene_symbol FROM genes)

            UNION ALL

            SELECT DISTINCT
                tissue AS Id,
                tissue AS Label,
                'tissue' AS node_type,
                0, NULL, NULL, NULL, NULL::INTEGER
            FROM expression

            UNION ALL

            SELECT DISTINCT
                'DISEASE:' || disease_name AS Id,
                disease_name AS Label,
                'disease' AS node_type,
                0, NULL, NULL, NULL, NULL::INTEGER
            FROM variants
            WHERE disease_name IS NOT NULL AND trim(disease_name) <> ''

            UNION ALL

            SELECT
                '{PANEL_DISEASE_ID}' AS Id,
                '{PANEL_DISEASE_LABEL}' AS Label,
                'disease' AS node_type,
                0, NULL, NULL,
                'ALS seed gene panel (reference disease)',
                NULL::INTEGER
            WHERE NOT EXISTS (
                SELECT 1 FROM variants
                WHERE disease_name = '{PANEL_DISEASE_LABEL}'
            )

            UNION ALL

            SELECT DISTINCT
                id AS Id,
                COALESCE(name, id) AS Label,
                'annotation' AS node_type,
                0, NULL, NULL, NULL, NULL::INTEGER
            FROM pathways_and_domains

            {foldseek_nodes}
        )
        SELECT
            Id,
            max(Label) AS Label,
            max(node_type) AS node_type,
            max(als_seed) AS als_seed,
            max(ensembl_id) AS ensembl_id,
            max(uniprot_id) AS uniprot_id,
            max(protein_description) AS protein_description,
            max(clinvar_row_count) AS clinvar_row_count
        FROM combined
        GROUP BY Id
    """

    e_path = os.path.join(out_dir, "atlas_unified_edges.csv")
    n_path = os.path.join(out_dir, "atlas_unified_nodes.csv")
    ne = _export_csv(conn, edges_sql, e_path)
    nn = _export_csv(conn, nodes_sql, n_path)
    print(f"unified: {ne} edges -> {e_path}")
    print(f"unified: {nn} nodes -> {n_path}")
    print("  Gephi: import atlas_unified_edges.csv then atlas_unified_nodes.csv")
    print("  Filter edges: Filters -> Attributes -> Equal -> network -> pick layers")
    print("  Highlight seeds: Appearance -> Nodes -> Partition -> als_seed")
    print("  Edge Weight is min-max normalized to [0,1] per network; see WeightRaw for originals")


SPLIT_EXPORTERS = {
    "string": lambda c, o, args: export_string(c, o),
    "gene-pathway": lambda c, o, args: export_gene_pathway(c, o),
    "foldseek": lambda c, o, args: export_foldseek(c, o, args.min_probability),
    "gene-bodyimpact": lambda c, o, args: export_gene_bodyimpact(c, o),
    "gene-disease": lambda c, o, args: export_gene_disease(c, o),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="Export ALS atlas DuckDB graphs for Gephi")
    parser.add_argument("--db", default=DEFAULT_DB, help="Path to als_genomic_atlas.duckdb")
    parser.add_argument("--out-dir", default=DEFAULT_OUT, help="Output directory for CSV files")
    parser.add_argument(
        "--network",
        choices=list(SPLIT_EXPORTERS) + ["unified", "all"],
        default="unified",
        help="unified = one combined graph (default); all = unified + per-layer CSVs",
    )
    parser.add_argument(
        "--split",
        action="store_true",
        help="With --network unified, also write per-layer CSV files",
    )
    parser.add_argument(
        "--min-probability",
        type=float,
        default=None,
        help="Foldseek: drop edges with probability below this (0–1)",
    )
    parser.add_argument(
        "--no-foldseek",
        action="store_true",
        help="Unified graph only: omit foldseek layer (smaller, easier to layout)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.db):
        print(f"Database not found: {args.db}", file=sys.stderr)
        return 1

    conn = duckdb.connect(args.db, read_only=True)
    try:
        if args.network in ("unified", "all"):
            export_unified(
                conn,
                args.out_dir,
                args.min_probability,
                include_foldseek=not args.no_foldseek,
            )
        if args.network == "all" or args.split:
            for name, fn in SPLIT_EXPORTERS.items():
                if args.no_foldseek and name == "foldseek":
                    continue
                fn(conn, args.out_dir, args)
        elif args.network in SPLIT_EXPORTERS:
            SPLIT_EXPORTERS[args.network](conn, args.out_dir, args)
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())