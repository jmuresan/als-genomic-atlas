import duckdb

def create_tables(conn: duckdb.DuckDBPyConnection):
    # Ingestion log table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS ingestion_log (
        source_name VARCHAR,
        query_params VARCHAR,
        status VARCHAR,
        record_count INTEGER,
        cache_path VARCHAR,
        error_message VARCHAR,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """)

    # 1. Genes table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS genes (
        gene_symbol VARCHAR PRIMARY KEY,
        ensembl_id VARCHAR,
        uniprot_id VARCHAR,
        chromosome VARCHAR,
        start_pos INTEGER,
        end_pos INTEGER,
        protein_description VARCHAR
    )
    """)

    # Transcripts table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS transcripts (
        transcript_id VARCHAR PRIMARY KEY,
        gene_symbol VARCHAR,
        mane_select BOOLEAN,
        length INTEGER,
        exons INTEGER,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 2. Variants table (combining ClinVar, dbSNP, gnomAD and AlphaGenome predictions)
    conn.execute("""
    CREATE TABLE IF NOT EXISTS variants (
        variant_id VARCHAR PRIMARY KEY,
        gene_symbol VARCHAR,
        clinical_significance VARCHAR,
        disease_name VARCHAR,
        rsid VARCHAR,
        chromosome VARCHAR,
        position INTEGER,
        hgvs VARCHAR,
        gnomad_pli DOUBLE,
        gnomad_loeuf DOUBLE,
        gnomad_allele_freq DOUBLE,
        alphagenome_consequence VARCHAR,
        alphagenome_pathogenicity DOUBLE,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 3. Regulatory Elements table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS regulatory_elements (
        element_id VARCHAR PRIMARY KEY,
        gene_symbol VARCHAR,
        element_type VARCHAR, -- promoter or enhancer
        score DOUBLE,
        ucsc_conservation_score DOUBLE,
        tfbs VARCHAR, -- comma-separated list
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 4. Expression table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS expression (
        gene_symbol VARCHAR,
        tissue VARCHAR,
        tpm DOUBLE,
        hpa_localization VARCHAR,
        hpa_score VARCHAR,
        PRIMARY KEY (gene_symbol, tissue),
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 5. Pathways & Domains table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS pathways_and_domains (
        id VARCHAR, -- pathway_id or domain_id
        gene_symbol VARCHAR,
        type VARCHAR, -- pathway, go_term, domain
        name VARCHAR,
        PRIMARY KEY (id, gene_symbol),
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 6. Interactions table
    #
    # STRING returns per-edge evidence channels in addition to the combined
    # `score`. We persist the experimental (escore), database/curated (dscore)
    # and text-mining (tscore) sub-scores so downstream consumers can tell a
    # physically-supported edge from a co-mention artifact. The combined
    # `confidence_score` alone CONFLATES these — a high combined score driven
    # purely by text-mining (escore~0, tscore~1) is a literature co-mention,
    # not evidence of a physical/functional interaction.
    conn.execute("""
    CREATE TABLE IF NOT EXISTS interactions (
        gene_a VARCHAR,
        gene_b VARCHAR,
        confidence_score DOUBLE,
        escore DOUBLE,   -- STRING experimental channel (physical evidence)
        dscore DOUBLE,   -- STRING database/curated channel
        tscore DOUBLE,   -- STRING text-mining channel (co-mention; NOT physical)
        PRIMARY KEY (gene_a, gene_b)
    )
    """)

    # 7. Clinical Trials & Drugs table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS clinical_trials_and_drugs (
        id VARCHAR, -- drug_id or trial_id
        gene_symbol VARCHAR,
        type VARCHAR, -- drug or trial
        name_or_title VARCHAR,
        max_clinical_phase DOUBLE,
        mechanism_of_action VARCHAR,
        status VARCHAR,
        PRIMARY KEY (id, gene_symbol),
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 8. Structures table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS structures (
        structure_id VARCHAR, -- pdb code or AF id
        gene_symbol VARCHAR,
        uniprot_id VARCHAR,
        type VARCHAR, -- PDB or AlphaFold
        plddt DOUBLE,
        disorder_score DOUBLE,
        method VARCHAR,
        PRIMARY KEY (structure_id, gene_symbol),
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 9. Foldseek matches table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS foldseek_matches (
        query_gene_symbol VARCHAR,
        target_id VARCHAR,
        db VARCHAR,
        probability DOUBLE,
        query_coverage DOUBLE,
        evalue DOUBLE,
        seq_identity DOUBLE,
        alignment_length INTEGER,
        PRIMARY KEY (query_gene_symbol, target_id),
        FOREIGN KEY (query_gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 10. Foldseek matched target drugs and trials table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS foldseek_matched_drugs_trials (
        query_gene_symbol VARCHAR,
        target_id VARCHAR,
        drug_or_trial_id VARCHAR,
        type VARCHAR, -- drug or trial
        name_or_title VARCHAR,
        max_clinical_phase DOUBLE,
        mechanism_of_action VARCHAR,
        status VARCHAR,
        purpose VARCHAR, -- what the drug/trial is for
        PRIMARY KEY (query_gene_symbol, target_id, drug_or_trial_id),
        FOREIGN KEY (query_gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 11. Foldseek similar compounds table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS foldseek_similar_compounds (
        query_gene_symbol VARCHAR,
        target_id VARCHAR,
        original_drug_id VARCHAR,
        similar_drug_id VARCHAR,
        name VARCHAR,
        similarity DOUBLE,
        max_clinical_phase DOUBLE,
        purpose VARCHAR,
        PRIMARY KEY (query_gene_symbol, target_id, original_drug_id, similar_drug_id),
        FOREIGN KEY (query_gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)


