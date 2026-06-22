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
        id VARCHAR PRIMARY KEY, -- pathway_id or domain_id
        gene_symbol VARCHAR,
        type VARCHAR, -- pathway, go_term, domain
        name VARCHAR,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 6. Interactions table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS interactions (
        gene_a VARCHAR,
        gene_b VARCHAR,
        confidence_score DOUBLE,
        PRIMARY KEY (gene_a, gene_b)
    )
    """)

    # 7. Clinical Trials & Drugs table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS clinical_trials_and_drugs (
        id VARCHAR PRIMARY KEY, -- drug_id or trial_id
        gene_symbol VARCHAR,
        type VARCHAR, -- drug or trial
        name_or_title VARCHAR,
        max_clinical_phase DOUBLE,
        mechanism_of_action VARCHAR,
        status VARCHAR,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)

    # 8. Structures table
    conn.execute("""
    CREATE TABLE IF NOT EXISTS structures (
        structure_id VARCHAR PRIMARY KEY, -- pdb code or AF id
        gene_symbol VARCHAR,
        uniprot_id VARCHAR,
        type VARCHAR, -- PDB or AlphaFold
        plddt DOUBLE,
        disorder_score DOUBLE,
        method VARCHAR,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    )
    """)
