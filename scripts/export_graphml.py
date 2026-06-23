#!/usr/bin/env python3
"""
Export all tables from DuckDB into a unified GraphML file.
This includes genes, transcripts, variants, regulatory elements, expression,
pathways/domains, interactions, clinical trials, drugs, structures, and foldseek matches.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import duckdb
import networkx as nx

DEFAULT_DB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "processed",
    "als_genomic_atlas.duckdb",
)
DEFAULT_OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "outputs",
    "als_genomic_atlas.graphml",
)

def clean_attrs(attrs: dict) -> dict:
    """Filter out None, empty strings, and convert complex types to strings."""
    cleaned = {}
    for k, v in attrs.items():
        if v is None or v == "":
            continue
        if isinstance(v, (int, float, bool, str)):
            cleaned[k] = v
        else:
            cleaned[k] = str(v)
    return cleaned

def add_node_safely(G: nx.DiGraph, node_id: str, node_type: str, label: str, attrs: dict):
    """Add a node or merge its attributes if it already exists."""
    full_attrs = {"node_type": node_type, "label": label, **attrs}
    cleaned = clean_attrs(full_attrs)
    if G.has_node(node_id):
        # Merge, preferring non-null/non-empty values
        for k, v in cleaned.items():
            G.nodes[node_id][k] = v
    else:
        G.add_node(node_id, **cleaned)

def add_edge_safely(G: nx.DiGraph, source: str, target: str, edge_type: str, network: str, attrs: dict):
    """Add an edge or merge its attributes if it already exists."""
    # Ensure source and target nodes exist (even as generic nodes if not added yet)
    if not G.has_node(source):
        G.add_node(source, node_type="unknown", label=source)
    if not G.has_node(target):
        G.add_node(target, node_type="unknown", label=target)
        
    full_attrs = {"edge_type": edge_type, "network": network, **attrs}
    cleaned = clean_attrs(full_attrs)
    if G.has_edge(source, target):
        for k, v in cleaned.items():
            G.edges[source, target][k] = v
    else:
        G.add_edge(source, target, **cleaned)

def export_graphml(db_path: str, out_path: str) -> nx.DiGraph:
    """Extracts all tables from the specified DuckDB database and saves as GraphML."""
    conn = duckdb.connect(db_path, read_only=True)
    G = nx.DiGraph()

    try:
        # 1. Genes (Seed Panel)
        genes = conn.execute(
            "SELECT gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description FROM genes"
        ).fetchall()
        for gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description in genes:
            add_node_safely(
                G,
                node_id=gene_symbol,
                node_type="gene",
                label=gene_symbol,
                attrs={
                    "ensembl_id": ensembl_id,
                    "uniprot_id": uniprot_id,
                    "chromosome": chromosome,
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                    "protein_description": protein_description,
                    "als_seed": 1,
                }
            )

        # 2. Transcripts
        transcripts = conn.execute(
            "SELECT transcript_id, gene_symbol, mane_select, length, exons FROM transcripts"
        ).fetchall()
        for transcript_id, gene_symbol, mane_select, length, exons in transcripts:
            add_node_safely(
                G,
                node_id=transcript_id,
                node_type="transcript",
                label=transcript_id,
                attrs={
                    "mane_select": bool(mane_select) if mane_select is not None else None,
                    "length": length,
                    "exons": exons,
                }
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=transcript_id,
                edge_type="has_transcript",
                network="gene-transcript",
                attrs={}
            )

        # 3. Variants & Diseases
        variants = conn.execute("""
            SELECT variant_id, gene_symbol, clinical_significance, disease_name, rsid, 
                   chromosome, position, hgvs, gnomad_pli, gnomad_loeuf, gnomad_allele_freq, 
                   alphagenome_consequence, alphagenome_pathogenicity 
            FROM variants
        """).fetchall()
        for (variant_id, gene_symbol, clinical_significance, disease_name, rsid, 
             chromosome, position, hgvs, gnomad_pli, gnomad_loeuf, gnomad_allele_freq, 
             alphagenome_consequence, alphagenome_pathogenicity) in variants:
            
            add_node_safely(
                G,
                node_id=variant_id,
                node_type="variant",
                label=variant_id,
                attrs={
                    "clinical_significance": clinical_significance,
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "hgvs": hgvs,
                    "gnomad_pli": gnomad_pli,
                    "gnomad_loeuf": gnomad_loeuf,
                    "gnomad_allele_freq": gnomad_allele_freq,
                    "alphagenome_consequence": alphagenome_consequence,
                    "alphagenome_pathogenicity": alphagenome_pathogenicity,
                }
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=variant_id,
                edge_type="has_variant",
                network="gene-variant",
                attrs={}
            )

            if disease_name and disease_name.strip():
                disease_id = f"DISEASE:{disease_name.strip()}"
                add_node_safely(
                    G,
                    node_id=disease_id,
                    node_type="disease",
                    label=disease_name.strip(),
                    attrs={}
                )
                add_edge_safely(
                    G,
                    source=variant_id,
                    target=disease_id,
                    edge_type="associated_with_disease",
                    network="variant-disease",
                    attrs={}
                )

        # Direct gene-disease associations
        gene_disease_clinvar = conn.execute("""
            SELECT gene_symbol, 'DISEASE:' || disease_name AS disease_id, disease_name, COUNT(*)::DOUBLE AS weight
            FROM variants
            WHERE disease_name IS NOT NULL AND trim(disease_name) <> ''
            GROUP BY gene_symbol, disease_name
        """).fetchall()
        for gene_symbol, disease_id, disease_name, weight in gene_disease_clinvar:
            add_node_safely(
                G,
                node_id=disease_id,
                node_type="disease",
                label=disease_name,
                attrs={}
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=disease_id,
                edge_type="associated_with_disease_direct",
                network="gene-disease",
                attrs={"Weight": weight, "source_type": "clinvar"}
            )

        # ALS phenotypic series panel connections
        panel_disease_id = "DISEASE:Amyotrophic lateral sclerosis"
        panel_disease_label = "Amyotrophic lateral sclerosis"
        add_node_safely(
            G,
            node_id=panel_disease_id,
            node_type="disease",
            label=panel_disease_label,
            attrs={"protein_description": "ALS seed gene panel (reference disease)"}
        )
        genes_no_variants = conn.execute("""
            SELECT gene_symbol FROM genes g
            WHERE NOT EXISTS (SELECT 1 FROM variants v WHERE v.gene_symbol = g.gene_symbol)
        """).fetchall()
        for (gene_symbol,) in genes_no_variants:
            add_edge_safely(
                G,
                source=gene_symbol,
                target=panel_disease_id,
                edge_type="seed_panel_disease",
                network="gene-disease",
                attrs={"Weight": 1.0, "source_type": "seed_panel"}
            )

        # 4. Regulatory Elements
        regulatory = conn.execute(
            "SELECT element_id, gene_symbol, element_type, score, ucsc_conservation_score, tfbs FROM regulatory_elements"
        ).fetchall()
        for element_id, gene_symbol, element_type, score, ucsc_conservation_score, tfbs in regulatory:
            add_node_safely(
                G,
                node_id=element_id,
                node_type="regulatory_element",
                label=element_id,
                attrs={
                    "element_type": element_type,
                    "score": score,
                    "ucsc_conservation_score": ucsc_conservation_score,
                    "tfbs": tfbs,
                }
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=element_id,
                edge_type="regulated_by",
                network="gene-regulatory",
                attrs={}
            )

        # 5. Expression
        expression = conn.execute(
            "SELECT gene_symbol, tissue, tpm, hpa_localization, hpa_score FROM expression"
        ).fetchall()
        for gene_symbol, tissue, tpm, hpa_localization, hpa_score in expression:
            add_node_safely(
                G,
                node_id=tissue,
                node_type="tissue",
                label=tissue,
                attrs={}
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=tissue,
                edge_type="expressed_in",
                network="gene-bodyimpact",
                attrs={
                    "tpm": tpm,
                    "hpa_localization": hpa_localization,
                    "hpa_score": hpa_score,
                }
            )

        # 6. Pathways and Domains (Annotations)
        pathways = conn.execute(
            "SELECT id, gene_symbol, type, name FROM pathways_and_domains"
        ).fetchall()
        for term_id, gene_symbol, term_type, name in pathways:
            add_node_safely(
                G,
                node_id=term_id,
                node_type="annotation",
                label=name if name else term_id,
                attrs={"annotation_type": term_type}
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=term_id,
                edge_type="associated_with",
                network="gene-pathway",
                attrs={}
            )

        # 7. STRING Interactions (PPI)
        interactions = conn.execute(
            "SELECT gene_a, gene_b, confidence_score FROM interactions"
        ).fetchall()
        for gene_a, gene_b, confidence_score in interactions:
            for g in [gene_a, gene_b]:
                if not G.has_node(g):
                    add_node_safely(G, node_id=g, node_type="gene", label=g, attrs={"als_seed": 0})
            
            src, tgt = sorted([gene_a, gene_b])
            add_edge_safely(
                G,
                source=src,
                target=tgt,
                edge_type="interacts_with",
                network="string",
                attrs={
                    "confidence_score": confidence_score,
                    "Weight": confidence_score,
                }
            )

        # 8. Clinical Trials and Drugs
        trials_drugs = conn.execute(
            "SELECT id, gene_symbol, type, name_or_title, max_clinical_phase, mechanism_of_action, status FROM clinical_trials_and_drugs"
        ).fetchall()
        for item_id, gene_symbol, item_type, name_or_title, max_clinical_phase, mechanism_of_action, status in trials_drugs:
            add_node_safely(
                G,
                node_id=item_id,
                node_type=item_type,
                label=name_or_title if name_or_title else item_id,
                attrs={
                    "max_clinical_phase": max_clinical_phase,
                    "mechanism_of_action": mechanism_of_action,
                    "status": status,
                }
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=item_id,
                edge_type="targeted_by" if item_type == "drug" else "involved_in",
                network="gene-drug-trial",
                attrs={
                    "max_clinical_phase": max_clinical_phase,
                    "mechanism_of_action": mechanism_of_action,
                    "status": status,
                }
            )

        # 9. Structures
        structures = conn.execute(
            "SELECT structure_id, gene_symbol, uniprot_id, type, plddt, disorder_score, method FROM structures"
        ).fetchall()
        for structure_id, gene_symbol, uniprot_id, structure_type, plddt, disorder_score, method in structures:
            add_node_safely(
                G,
                node_id=structure_id,
                node_type="structure",
                label=structure_id,
                attrs={
                    "structure_type": structure_type,
                    "uniprot_id": uniprot_id,
                    "plddt": plddt,
                    "disorder_score": disorder_score,
                    "method": method,
                }
            )
            add_edge_safely(
                G,
                source=gene_symbol,
                target=structure_id,
                edge_type="has_structure",
                network="gene-structure",
                attrs={}
            )

        # 10. Foldseek Matches
        foldseek = conn.execute(
            "SELECT query_gene_symbol, target_id, db, probability, query_coverage, evalue, seq_identity, alignment_length FROM foldseek_matches"
        ).fetchall()
        for query_gene_symbol, target_id, db, probability, query_coverage, evalue, seq_identity, alignment_length in foldseek:
            match = re.search(r'(AF-[A-Z0-9]+)', target_id)
            label = match.group(1) if match else target_id[:80]
            
            add_node_safely(
                G,
                node_id=target_id,
                node_type="foldseek_target",
                label=label,
                attrs={"db": db}
            )
            add_edge_safely(
                G,
                source=query_gene_symbol,
                target=target_id,
                edge_type="structurally_similar_to",
                network="foldseek",
                attrs={
                    "db": db,
                    "probability": probability,
                    "query_coverage": query_coverage,
                    "evalue": evalue,
                    "seq_identity": seq_identity,
                    "alignment_length": alignment_length,
                }
            )

        # 11. Foldseek Matched Drugs and Trials
        fs_drugs = conn.execute(
            "SELECT query_gene_symbol, target_id, drug_or_trial_id, type, name_or_title, max_clinical_phase, mechanism_of_action, status, purpose FROM foldseek_matched_drugs_trials"
        ).fetchall()
        for query_gene_symbol, target_id, drug_or_trial_id, item_type, name_or_title, max_clinical_phase, mechanism_of_action, status, purpose in fs_drugs:
            add_node_safely(
                G,
                node_id=drug_or_trial_id,
                node_type=item_type,
                label=name_or_title if name_or_title else drug_or_trial_id,
                attrs={
                    "max_clinical_phase": max_clinical_phase,
                    "mechanism_of_action": mechanism_of_action,
                    "status": status,
                    "purpose": purpose,
                }
            )
            add_edge_safely(
                G,
                source=target_id,
                target=drug_or_trial_id,
                edge_type="targeted_by" if item_type == "drug" else "involved_in",
                network="foldseek-drug-trial",
                attrs={
                    "max_clinical_phase": max_clinical_phase,
                    "mechanism_of_action": mechanism_of_action,
                    "status": status,
                    "purpose": purpose,
                }
            )

        # 12. Foldseek Similar Compounds
        similar_compounds = conn.execute(
            "SELECT query_gene_symbol, target_id, original_drug_id, similar_drug_id, name, similarity, max_clinical_phase, purpose FROM foldseek_similar_compounds"
        ).fetchall()
        for query_gene_symbol, target_id, original_drug_id, similar_drug_id, name, similarity, max_clinical_phase, purpose in similar_compounds:
            if not G.has_node(original_drug_id):
                add_node_safely(G, node_id=original_drug_id, node_type="drug", label=original_drug_id, attrs={})
                
            add_node_safely(
                G,
                node_id=similar_drug_id,
                node_type="drug",
                label=name if name else similar_drug_id,
                attrs={
                    "max_clinical_phase": max_clinical_phase,
                    "purpose": purpose,
                }
            )
            add_edge_safely(
                G,
                source=original_drug_id,
                target=similar_drug_id,
                edge_type="similar_to",
                network="drug-similarity",
                attrs={
                    "similarity": similarity,
                    "max_clinical_phase": max_clinical_phase,
                    "purpose": purpose,
                }
            )

        # Save to GraphML
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        nx.write_graphml(G, out_path)

    finally:
        conn.close()

    return G

def main() -> int:
    parser = argparse.ArgumentParser(description="Export ALS atlas DuckDB database to GraphML")
    parser.add_argument("--db", default=DEFAULT_DB, help="Path to als_genomic_atlas.duckdb")
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output path for GraphML file")
    args = parser.parse_args()

    if not os.path.isfile(args.db):
        print(f"Database not found: {args.db}", file=sys.stderr)
        return 1

    print(f"Connecting to database: {args.db}")
    G = export_graphml(args.db, args.out)

    # Print statistics
    print("\n--- Export Statistics ---")
    print(f"Total Nodes: {G.number_of_nodes()}")
    print(f"Total Edges: {G.number_of_edges()}")
    
    node_types = {}
    for n, data in G.nodes(data=True):
        nt = data.get("node_type", "unknown")
        node_types[nt] = node_types.get(nt, 0) + 1
    print("\nNodes by type:")
    for nt, count in sorted(node_types.items()):
        print(f"  {nt}: {count}")

    edge_types = {}
    for u, v, data in G.edges(data=True):
        et = data.get("edge_type", "unknown")
        edge_types[et] = edge_types.get(et, 0) + 1
    print("\nEdges by type:")
    for et, count in sorted(edge_types.items()):
        print(f"  {et}: {count}")

    print(f"\nWriting GraphML file to: {args.out}")
    print("Export completed successfully!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
