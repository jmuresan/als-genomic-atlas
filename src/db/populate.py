import json
import duckdb
import logging
from typing import Dict, Any, Optional, List

logger = logging.getLogger("als_atlas.populate")

STAGE_TO_PHASE = {
    "APPROVAL": 4.0, "PHASE_IV": 4.0, "PHASE_III": 3.0, "PHASE_II": 2.0,
    "PHASE_I": 1.0, "EARLY_PHASE_I": 0.5, "PRECLINICAL": 0.0, "PHASE_0": 0.0,
}

def stage_to_phase(stage: Optional[str]) -> Optional[float]:
    if not stage:
        return None
    return STAGE_TO_PHASE.get(str(stage).strip().upper())

def log_ingestion(conn: duckdb.DuckDBPyConnection, source_name: str, query_params: Optional[Dict[str, Any]], status: str, record_count: int, cache_path: Optional[str], error_message: Optional[str]):
    serialized_params = json.dumps(query_params) if query_params else "{}"
    if error_message:
        error_message = error_message[:1000]
    conn.execute("""
    INSERT INTO ingestion_log (source_name, query_params, status, record_count, cache_path, error_message)
    VALUES (?, ?, ?, ?, ?, ?)
    """, [source_name, serialized_params, status, record_count, cache_path, error_message])

def populate_gene(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    uniprot_id = data.get("uniprot_id")
    
    ensembl = data.get("ensembl", {}) or {}
    ensembl_id = ensembl.get("ensembl_id")
    coords = ensembl.get("coordinates", {}) or {}
    chromosome = coords.get("chr")
    start_pos = coords.get("start")
    end_pos = coords.get("end")
    
    uniprot = data.get("uniprot", {}) or {}
    protein_description = None
    desc = uniprot.get("proteinDescription", {}) or {}
    if "recommendedName" in desc:
        protein_description = desc["recommendedName"].get("fullName", {}).get("value")
        
    conn.execute("""
    INSERT OR REPLACE INTO genes (gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    """, [gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description])

    # Transcripts
    transcripts = ensembl.get("transcripts", []) or []
    for tx in transcripts:
        tx_id = tx.get("id")
        mane = tx.get("mane_select", False)
        length = tx.get("length")
        exons = tx.get("exons")
        if tx_id:
            conn.execute("""
            INSERT OR REPLACE INTO transcripts (transcript_id, gene_symbol, mane_select, length, exons)
            VALUES (?, ?, ?, ?, ?)
            """, [tx_id, gene_symbol, mane, length, exons])

def _clinvar_grch38_loc(info: Dict[str, Any]):
    """Per-variant GRCh38 (chromosome, position) from a ClinVar esummary record.

    ClinVar records carry their own authoritative coordinates in
    variation_set[].variation_loc[]; prefer the GRCh38 assembly entry. Returns
    (None, None) if unavailable so the caller can fall back.
    """
    for vset in (info.get("variation_set") or []):
        for loc in (vset.get("variation_loc") or []):
            if loc.get("assembly_name") == "GRCh38" and loc.get("chr"):
                start = loc.get("start") or loc.get("display_start")
                try:
                    pos = int(start) if start not in (None, "") else None
                except (TypeError, ValueError):
                    pos = None
                return loc.get("chr"), pos
    return None, None


def _clinvar_rsid(info: Dict[str, Any]):
    """Per-variant dbSNP rsid from a ClinVar record's variation_xrefs, if present."""
    for vset in (info.get("variation_set") or []):
        for xref in (vset.get("variation_xrefs") or []):
            if (xref.get("db_source") or "").lower() == "dbsnp":
                rid = xref.get("db_id")
                if rid:
                    rid = str(rid)
                    return rid if rid.startswith("rs") else "rs" + rid
    return None


def populate_variants(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")

    # ClinVar
    clinvar = data.get("clinvar", {}) or {}
    clinvar_result = clinvar.get("result", {}) or {}

    # dbSNP — gene-level fallback ONLY. Historically these per-gene values were
    # written onto EVERY ClinVar variant of the gene, collapsing all coordinates
    # to one (corrupt) locus. Real coordinates are now taken per-variant from the
    # ClinVar record below; dbSNP is used only when ClinVar lacks a coordinate.
    dbsnp = data.get("dbsnp", {}) or {}
    fallback_rsid = dbsnp.get("rsid")
    fallback_chromosome = dbsnp.get("chromosome")
    fallback_position = dbsnp.get("position")
    fallback_hgvs = dbsnp.get("hgvs")
    
    # gnomAD
    gnomad = data.get("gnomad", {}) or {}
    pli = gnomad.get("pli")
    loeuf = gnomad.get("loeuf")
    allele_freq = gnomad.get("allele_freq")
    
    # AlphaGenome
    alphagenome = data.get("alphagenome", {}) or {}
    consequence = alphagenome.get("non_coding_variant_consequence")
    pathogenicity = alphagenome.get("pathogenicity_score")
    
    variant_found = False
    for uid, info in clinvar_result.items():
        if uid == "uids":
            continue
        variant_id = info.get("accession") or info.get("uid")
        if not variant_id:
            continue
        
        germline = info.get("germline_classification") or {}
        clinical_sig = germline.get("description")
        traits = germline.get("trait_set") or []
        disease_name = traits[0].get("trait_name") if traits else None

        # Per-variant coordinates from the ClinVar record (authoritative),
        # falling back to the gene-level dbSNP record only if absent.
        v_chr, v_pos = _clinvar_grch38_loc(info)
        v_chr = v_chr if v_chr is not None else fallback_chromosome
        v_pos = v_pos if v_pos is not None else fallback_position
        v_hgvs = info.get("title") or fallback_hgvs
        v_rsid = _clinvar_rsid(info) or fallback_rsid

        conn.execute("""
        INSERT OR REPLACE INTO variants (
            variant_id, gene_symbol, clinical_significance, disease_name,
            rsid, chromosome, position, hgvs, gnomad_pli, gnomad_loeuf,
            gnomad_allele_freq, alphagenome_consequence, alphagenome_pathogenicity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            variant_id, gene_symbol, clinical_sig, disease_name,
            v_rsid, v_chr, v_pos, v_hgvs, pli, loeuf,
            allele_freq, consequence, pathogenicity
        ])
        variant_found = True

    # De-mocked: the original inserted a fabricated "VAR_<gene>_MOCK"
    # Benign/Likely Benign row whenever a gene had no real ClinVar hit. That
    # invented data, so it has been removed. Genes with no ClinVar variants
    # simply contribute no variant rows (their gnomAD/dbSNP signal still lands
    # via genes/other tables).
    _ = variant_found

def populate_regulatory_elements(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    
    encode = data.get("encode", {}) or {}
    ucsc = data.get("ucsc", {}) or {}
    
    promoters = encode.get("promoters", []) or []
    enhancers = encode.get("enhancers", []) or []
    
    conservation = ucsc.get("conservation_score", 0.0)
    tfbs_list = ",".join(ucsc.get("tfbs", []))
    
    for p in promoters:
        conn.execute("""
        INSERT OR REPLACE INTO regulatory_elements (element_id, gene_symbol, element_type, score, ucsc_conservation_score, tfbs)
        VALUES (?, ?, 'promoter', ?, ?, ?)
        """, [p.get("id"), gene_symbol, p.get("score"), conservation, tfbs_list])
        
    for e in enhancers:
        conn.execute("""
        INSERT OR REPLACE INTO regulatory_elements (element_id, gene_symbol, element_type, score, ucsc_conservation_score, tfbs)
        VALUES (?, ?, 'enhancer', ?, ?, ?)
        """, [e.get("id"), gene_symbol, e.get("score"), conservation, tfbs_list])

def populate_expression(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    gtex = data.get("gtex", {}) or {}
    hpa = data.get("hpa", {}) or {}
    
    tissues = gtex.get("tissues", []) or []
    localization = hpa.get("localization")
    hpa_score = hpa.get("score")
    
    for t in tissues:
        tissue_name = t.get("tissue")
        tpm = t.get("tpm")
        if tissue_name:
            conn.execute("""
            INSERT OR REPLACE INTO expression (gene_symbol, tissue, tpm, hpa_localization, hpa_score)
            VALUES (?, ?, ?, ?, ?)
            """, [gene_symbol, tissue_name, tpm, localization, hpa_score])

def populate_pathways_and_domains(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    
    reactome = data.get("reactome", []) or []
    quickgo = data.get("quickgo", {}) or {}
    interpro = data.get("interpro", {}) or {}
    
    # Reactome pathways
    for p in reactome:
        pathway_id = p.get("stId")
        pathway_name = p.get("displayName")
        if pathway_id:
            conn.execute("""
            INSERT OR REPLACE INTO pathways_and_domains (id, gene_symbol, type, name)
            VALUES (?, ?, 'pathway', ?)
            """, [pathway_id, gene_symbol, pathway_name])
            
    # GO terms from QuickGO
    go_results = quickgo.get("results", []) or []
    for go in go_results:
        go_id = go.get("goId")
        go_name = go.get("goName")
        if go_id:
            conn.execute("""
            INSERT OR REPLACE INTO pathways_and_domains (id, gene_symbol, type, name)
            VALUES (?, ?, 'go_term', ?)
            """, [go_id, gene_symbol, go_name])
            
    # InterPro domains
    ip_results = interpro.get("results", []) or []
    for ip in ip_results:
        ip_id = ip.get("metadata", {}).get("accession")
        ip_name = ip.get("metadata", {}).get("name")
        if ip_id:
            conn.execute("""
            INSERT OR REPLACE INTO pathways_and_domains (id, gene_symbol, type, name)
            VALUES (?, ?, 'domain', ?)
            """, [ip_id, gene_symbol, ip_name])

def populate_interactions(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    string_list = data.get("string", []) or []
    for item in string_list:
        gene_a = item.get("preferredName_A")
        gene_b = item.get("preferredName_B")
        score = item.get("score")
        if gene_a and gene_b:
            gene_a, gene_b = sorted([gene_a, gene_b])
            # Persist the STRING evidence channels alongside the combined score
            # so consumers can distinguish physical/curated edges (escore/dscore)
            # from text-mining co-mentions (tscore). All come from the same
            # cached STRING response, so this stays deterministic on rebuild.
            conn.execute("""
            INSERT OR REPLACE INTO interactions (gene_a, gene_b, confidence_score, escore, dscore, tscore)
            VALUES (?, ?, ?, ?, ?, ?)
            """, [gene_a, gene_b, score, item.get("escore"), item.get("dscore"), item.get("tscore")])

def populate_clinical_trials_and_drugs(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    
    # Open Targets (disease associations and drugs)
    ot = data.get("open_targets", {}) or {}
    drugs = ot.get("drugs", []) or []
    for row in drugs:
        drug = row.get("drug", {}) or {}
        drug_id = drug.get("id")
        drug_name = drug.get("name")
        stage = drug.get("maximumClinicalStage") or row.get("maxClinicalStage")
        phase = stage_to_phase(stage)
        moa_rows = (drug.get("mechanismsOfAction") or {}).get("rows") or []
        mech = moa_rows[0].get("mechanismOfAction") if moa_rows else None
        
        if drug_id:
            conn.execute("""
            INSERT OR REPLACE INTO clinical_trials_and_drugs (id, gene_symbol, type, name_or_title, max_clinical_phase, mechanism_of_action)
            VALUES (?, ?, 'drug', ?, ?, ?)
            """, [drug_id, gene_symbol, drug_name, phase, mech])

    # Clinical Trials from ClinicalTrials.gov
    ct = data.get("clinical_trials", {}) or {}
    trials = ct.get("trials", []) or []
    for trial in trials:
        nct_id = trial.get("nct_id")
        title = trial.get("title")
        status = trial.get("status")
        if nct_id:
            conn.execute("""
            INSERT OR REPLACE INTO clinical_trials_and_drugs (id, gene_symbol, type, name_or_title, status)
            VALUES (?, ?, 'trial', ?, ?)
            """, [nct_id, gene_symbol, title, status])

def populate_structures(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    uniprot_id = data.get("uniprot_id")
    
    # AlphaFold
    af = data.get("alphafold", {}) or {}
    if isinstance(af, list):
        af = af[0] if af else {}
        
    plddt = af.get("globalMetricValue") or af.get("plddt")
    disorder = af.get("disorder_score")
    if plddt is not None:
        conn.execute("""
        INSERT OR REPLACE INTO structures (structure_id, gene_symbol, uniprot_id, type, plddt, disorder_score)
        VALUES (?, ?, ?, 'AlphaFold', ?, ?)
        """, [f"AF-{uniprot_id}", gene_symbol, uniprot_id, plddt, disorder])
        
    # PDB
    pdb = data.get("pdb", {}) or {}
    pdb_ids = pdb.get("pdb_ids", []) or []
    method = pdb.get("method")
    for pdb_id in pdb_ids:
        conn.execute("""
        INSERT OR REPLACE INTO structures (structure_id, gene_symbol, uniprot_id, type, method)
        VALUES (?, ?, ?, 'PDB', ?)
        """, [pdb_id, gene_symbol, uniprot_id, method])

def populate_foldseek_matches(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    matches = data.get("foldseek_matches", []) or []
    for m in matches:
        conn.execute("""
        INSERT OR REPLACE INTO foldseek_matches (
            query_gene_symbol, target_id, db, probability, query_coverage, evalue, seq_identity, alignment_length
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            gene_symbol, m.get("target_id"), m.get("db"), m.get("probability"),
            m.get("query_coverage"), m.get("eval"), m.get("seqId"), m.get("alnLength")
        ])

def populate_foldseek_matched_drugs_trials(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    drugs = data.get("foldseek_drugs", []) or []
    for d in drugs:
        conn.execute("""
        INSERT OR REPLACE INTO foldseek_matched_drugs_trials (
            query_gene_symbol, target_id, drug_or_trial_id, type, name_or_title, 
            max_clinical_phase, mechanism_of_action, status, purpose
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            gene_symbol, d.get("target_id"), d.get("drug_id"), d.get("type"), d.get("name_or_title"),
            d.get("max_clinical_phase"), d.get("mechanism_of_action"), d.get("status"), d.get("purpose")
        ])

def populate_foldseek_similar_compounds(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    gene_symbol = data.get("gene")
    compounds = data.get("similar_compounds", []) or []
    for c in compounds:
        conn.execute("""
        INSERT OR REPLACE INTO foldseek_similar_compounds (
            query_gene_symbol, target_id, original_drug_id, similar_drug_id, name,
            similarity, max_clinical_phase, purpose
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            gene_symbol, c.get("target_id"), c.get("original_drug_id"), c.get("similar_drug_id"),
            c.get("name"), c.get("similarity"), c.get("max_clinical_phase"), c.get("purpose")
        ])

def populate_all(conn: duckdb.DuckDBPyConnection, data: Dict[str, Any]):
    """Populates all 11 categories of biological information into DuckDB."""
    populate_gene(conn, data)
    populate_variants(conn, data)
    populate_regulatory_elements(conn, data)
    populate_expression(conn, data)
    populate_pathways_and_domains(conn, data)
    populate_interactions(conn, data)
    populate_clinical_trials_and_drugs(conn, data)
    populate_structures(conn, data)
    populate_foldseek_matches(conn, data)
    populate_foldseek_matched_drugs_trials(conn, data)
    populate_foldseek_similar_compounds(conn, data)
