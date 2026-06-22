import os
import duckdb
import logging
from typing import Dict, Any, List

logger = logging.getLogger("als_atlas.generate_atlas")

def generate_report(db_path: str, output_path: str):
    """Compiles the DuckDB database contents into a comprehensive Markdown Genomic Atlas."""
    logger.info(f"Generating Genomic Atlas report from {db_path}...")
    conn = duckdb.connect(db_path, read_only=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Fetch all genes
    genes = conn.execute("""
        SELECT gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description 
        FROM genes 
        ORDER BY gene_symbol
    """).fetchall()
    
    with open(output_path, "w", encoding="utf-8") as f:
        # Header
        f.write("# ALS Genomic Atlas\n\n")
        f.write("A publication-quality molecular and genomic landscape map of Amyotrophic Lateral Sclerosis (ALS) associated genes. Compiled dynamically from DuckDB containing multi-omics data integrated across 8 distinct categories.\n\n")
        
        f.write("## Table of Contents\n")
        for row in genes:
            symbol = row[0]
            f.write(f"- [{symbol}](#{symbol.lower()})\n")
        f.write("\n---\n\n")
        
        for row in genes:
            symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_desc = row
            
            f.write(f"## {symbol}\n\n")
            f.write(f"**Protein Description:** {protein_desc or 'N/A'}\n\n")
            f.write(f"### Category 1: Gene & Transcript Mapping\n")
            f.write(f"- **Ensembl ID:** `{ensembl_id or 'N/A'}`\n")
            f.write(f"- **UniProt ID:** `{uniprot_id or 'N/A'}`\n")
            f.write(f"- **Location:** Chr{chromosome or 'N/A'}:{start_pos or 'N/A'}-{end_pos or 'N/A'}\n\n")
            
            # Transcripts table
            transcripts = conn.execute("""
                SELECT transcript_id, mane_select, length, exons 
                FROM transcripts 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            if transcripts:
                f.write("#### Transcripts\n")
                f.write("| Transcript ID | MANE Select | Length (bp) | Exons |\n")
                f.write("| --- | --- | --- | --- |\n")
                for tx in transcripts:
                    f.write(f"| {tx[0]} | {'Yes' if tx[1] else 'No'} | {tx[2]} | {tx[3]} |\n")
                f.write("\n")
            
            # Category 2: Variants & Pathogenicity
            variants = conn.execute("""
                SELECT variant_id, clinical_significance, disease_name, rsid, gnomad_pli, gnomad_loeuf, gnomad_allele_freq, alphagenome_consequence, alphagenome_pathogenicity 
                FROM variants 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            f.write("### Category 2: Genetic Variants & Pathogenicity\n")
            if variants:
                f.write("| Variant ID | Clinical Significance | gnomAD pLI | gnomAD LOEUF | AlphaGenome Pred |\n")
                f.write("| --- | --- | --- | --- | --- |\n")
                for v in variants:
                    v_id, sig, disease, rs, pli, loeuf, af, cons, path = v
                    f.write(f"| `{v_id}` (rsid: {rs or 'N/A'}) | {sig or 'N/A'} | {pli if pli is not None else 'N/A'} | {loeuf if loeuf is not None else 'N/A'} | {cons or 'N/A'} (Path Score: {path if path is not None else 'N/A'}) |\n")
                f.write("\n")
            else:
                f.write("*No variants recorded.*\n\n")
                
            # Category 3: Transcriptional Regulation & Epigenomics
            reg_elements = conn.execute("""
                SELECT element_id, element_type, score, ucsc_conservation_score, tfbs 
                FROM regulatory_elements 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            f.write("### Category 3: Transcriptional Regulation & Epigenomics\n")
            if reg_elements:
                f.write("| Element ID | Type | Score | UCSC Conservation | TFBS Motifs |\n")
                f.write("| --- | --- | --- | --- | --- |\n")
                for re_el in reg_elements:
                    f.write(f"| `{re_el[0]}` | {re_el[1].capitalize()} | {re_el[2]} | {re_el[3]} | {re_el[4] or 'N/A'} |\n")
                f.write("\n")
            else:
                f.write("*No regulatory elements annotated.*\n\n")

            # Category 4: Expression & Tissue Specificity
            expression = conn.execute("""
                SELECT tissue, tpm, hpa_localization, hpa_score 
                FROM expression 
                WHERE gene_symbol = ?
                ORDER BY tpm DESC
            """, [symbol]).fetchall()
            
            f.write("### Category 4: Expression & Tissue Specificity\n")
            if expression:
                f.write(f"- **HPA Cellular Localization:** {expression[0][2] or 'N/A'} (HPA Score: {expression[0][3] or 'N/A'})\n\n")
                f.write("#### Tissue TPM Expression (GTEx)\n")
                f.write("| Tissue | Median TPM |\n")
                f.write("| --- | --- |\n")
                for exp in expression:
                    f.write(f"| {exp[0]} | {exp[1]:.2f} |\n")
                f.write("\n")
            else:
                f.write("*No expression data mapped.*\n\n")
                
            # Category 5: Pathways & Functional Annotation
            pathways = conn.execute("""
                SELECT id, type, name 
                FROM pathways_and_domains 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            f.write("### Category 5: Pathways & Functional Annotation\n")
            if pathways:
                f.write("| ID | Type | Name |\n")
                f.write("| --- | --- | --- |\n")
                for pw in pathways:
                    f.write(f"| `{pw[0]}` | {pw[1]} | {pw[2]} |\n")
                f.write("\n")
            else:
                f.write("*No pathway/GO annotations mapped.*\n\n")
                
            # Category 6: Network Interactions
            interactions = conn.execute("""
                SELECT gene_a, gene_b, confidence_score 
                FROM interactions 
                WHERE gene_a = ? OR gene_b = ?
                ORDER BY confidence_score DESC
            """, [symbol, symbol]).fetchall()
            
            f.write("### Category 6: Network Interactions (STRING)\n")
            if interactions:
                f.write("| Interactor A | Interactor B | Confidence Score |\n")
                f.write("| --- | --- | --- |\n")
                for inter in interactions:
                    f.write(f"| {inter[0]} | {inter[1]} | {inter[2]:.3f} |\n")
                f.write("\n")
            else:
                f.write("*No interactors mapped.*\n\n")
                
            # Category 7: Clinical Translation & Druggability
            drugs_trials = conn.execute("""
                SELECT id, type, name_or_title, max_clinical_phase, mechanism_of_action, status 
                FROM clinical_trials_and_drugs 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            f.write("### Category 7: Clinical Translation & Druggability\n")
            if drugs_trials:
                f.write("| ID | Type | Name/Title | Phase/Status | Details |\n")
                f.write("| --- | --- | --- | --- | --- |\n")
                for dt in drugs_trials:
                    if dt[1] == 'drug':
                        f.write(f"| `{dt[0]}` | Drug | {dt[2]} | Phase {dt[3]} | MoA: {dt[4] or 'N/A'} |\n")
                    else:
                        f.write(f"| `{dt[0]}` | Trial | {dt[2]} | {dt[5]} | N/A |\n")
                f.write("\n")
            else:
                f.write("*No drug/trial associations.*\n\n")
                
            # Category 8: 3D Structural Biology
            structures = conn.execute("""
                SELECT structure_id, type, plddt, disorder_score, method 
                FROM structures 
                WHERE gene_symbol = ?
            """, [symbol]).fetchall()
            
            f.write("### Category 8: 3D Structural Biology\n")
            if structures:
                f.write("| Structure ID | Database | pLDDT / Method | Disorder Score |\n")
                f.write("| --- | --- | --- | --- |\n")
                for s in structures:
                    if s[1] == 'AlphaFold':
                        f.write(f"| `{s[0]}` | AlphaFold | pLDDT: {s[2]} | Disorder: {s[3]} |\n")
                    else:
                        f.write(f"| `{s[0]}` | PDB | Method: {s[4]} | N/A |\n")
                f.write("\n")
            else:
                f.write("*No 3D structure mappings.*\n\n")
                
            # Category 9: Foldseek Structural Similarity Matches.
            # Ordered by E-value ascending (most statistically significant
            # alignments first); missing E-values sort last.
            foldseek = conn.execute("""
                SELECT target_id, db, probability, query_coverage, evalue, seq_identity, alignment_length
                FROM foldseek_matches
                WHERE query_gene_symbol = ?
                ORDER BY evalue ASC NULLS LAST
            """, [symbol]).fetchall()
            
            f.write("### Category 9: Foldseek Structural Similarity Matches\n")
            if foldseek:
                f.write("| Target ID | Database | Probability | Query Coverage | E-value | Seq Identity | Aln Length |\n")
                f.write("| --- | --- | --- | --- | --- | --- | --- |\n")
                for fs in foldseek:
                    target_id, db, prob, q_cov, evalue, seq_id, aln_len = fs
                    q_cov_str = f"{q_cov * 100:.1f}%" if q_cov is not None else "N/A"
                    prob_str = f"{prob:.3f}" if prob is not None else "N/A"
                    seq_id_str = f"{seq_id * 100:.1f}%" if seq_id is not None else "N/A"
                    f.write(f"| `{target_id}` | {db} | {prob_str} | {q_cov_str} | {evalue} | {seq_id_str} | {aln_len} |\n")
                f.write("\n")
            else:
                f.write("*No Foldseek structural similarity matches found.*\n\n")
                
            # Category 10: Matched Target Drugs & Clinical Trials
            # Ordered by clinical phase (then status), highest on top, so the
            # most advanced drugs/trials for matched targets surface first.
            matched_drugs = conn.execute("""
                SELECT target_id, drug_or_trial_id, name_or_title, type, max_clinical_phase, mechanism_of_action, status, purpose
                FROM foldseek_matched_drugs_trials
                WHERE query_gene_symbol = ?
                ORDER BY max_clinical_phase DESC NULLS LAST, status DESC NULLS LAST
            """, [symbol]).fetchall()
            
            f.write("### Category 10: Matched Target Drugs & Clinical Trials\n")
            if matched_drugs:
                f.write("| Matched Target | Drug/Trial ID | Name/Title | Type | Phase/Status | Purpose / Description |\n")
                f.write("| --- | --- | --- | --- | --- | --- |\n")
                for md in matched_drugs:
                    target_id, drug_id, name, d_type, phase, moa, status, purpose = md
                    phase_str = f"Phase {phase}" if phase is not None else (status or "N/A")
                    f.write(f"| `{target_id}` | `{drug_id}` | {name or 'N/A'} | {d_type.capitalize()} | {phase_str} | {purpose or 'N/A'} |\n")
                f.write("\n")
            else:
                f.write("*No drugs or clinical trials associated with matched target proteins.*\n\n")
                
            # Category 11: Similar Compounds & Repurposing Candidates
            # Ordered by clinical stage (max phase), highest on top, so the most
            # advanced repurposing candidates surface first; similarity breaks ties.
            similar_compounds = conn.execute("""
                SELECT target_id, original_drug_id, similar_drug_id, name, similarity, max_clinical_phase, purpose
                FROM foldseek_similar_compounds
                WHERE query_gene_symbol = ?
                ORDER BY max_clinical_phase DESC NULLS LAST, similarity DESC
            """, [symbol]).fetchall()
            
            f.write("### Category 11: Similar Compounds & Repurposing Candidates\n")
            if similar_compounds:
                f.write("| Matched Target | Original Drug | Similar Drug ID | Similar Drug Name | Similarity % | Max Phase | Indications & Mechanism |\n")
                f.write("| --- | --- | --- | --- | --- | --- | --- |\n")
                for sc in similar_compounds:
                    target_id, orig_id, sim_id, sim_name, similarity, max_phase, purpose = sc
                    sim_pct = f"{similarity:.1f}%" if similarity is not None else "N/A"
                    f.write(f"| `{target_id}` | `{orig_id}` | `{sim_id}` | {sim_name or 'N/A'} | {sim_pct} | Phase {max_phase} | {purpose or 'N/A'} |\n")
                f.write("\n")
            else:
                f.write("*No similar compounds with qualifying max clinical phase found.*\n\n")
                
            f.write("\n---\n\n")
            
    conn.close()
    logger.info(f"Atlas report written successfully to {output_path}")
