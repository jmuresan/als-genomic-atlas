import os
import argparse
import logging
import duckdb
import sys
from src.config import Config
from src.ingest.cache import DiskCache
from src.ingest.client import IngestionManager
from src.db.schema import create_tables
from src.db.populate import populate_all, log_ingestion
from src.atlas.generate_atlas import generate_report

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("als_atlas.pipeline")

def run_pipeline(config_path: str = None, test_run: bool = False):
    logger.info("Initializing ALS Genomic Atlas pipeline...")
    
    # 1. Load config
    config = Config(config_path)
    logger.info(f"Loaded config. Offline mode: {config.offline_mode}")
    
    # 2. Setup paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    db_dir = os.path.join(base_dir, "data", "processed")
    os.makedirs(db_dir, exist_ok=True)
    db_path = os.path.join(db_dir, "als_genomic_atlas.duckdb")
    
    # Connect and setup schema
    conn = duckdb.connect(db_path)
    create_tables(conn)
    
    # 3. Setup cache & Ingestion Manager
    cache = DiskCache(config.cache_dir, config.offline_mode)
    manager = IngestionManager(cache)
    
    # Determine gene list (limited to 5 for test_run/speed)
    genes_to_process = config.seed_genes
    if test_run:
        logger.info("Test run active: limiting to first 5 seed genes.")
        genes_to_process = config.seed_genes[:5]
        
    logger.info(f"Processing genes: {', '.join(genes_to_process)}")
    
    # 4. Ingest and populate
    for gene in genes_to_process:
        try:
            # Fetch data
            data = manager.fetch_all_data(gene)
            
            # Populate tables
            populate_all(conn, data)
            
            # Log success
            log_ingestion(conn, "ingestion_manager", {"gene": gene}, "SUCCESS", 1, config.cache_dir, None)
            logger.info(f"Successfully processed and populated data for {gene}")
        except Exception as e:
            logger.error(f"Failed to process gene {gene}: {e}")
            log_ingestion(conn, "ingestion_manager", {"gene": gene}, "FAILED", 0, None, str(e))
            if not test_run:
                # In standard pipeline, fail hard on fatal issues
                conn.close()
                raise e
                
    conn.close()
    
    # 5. Generate Atlas Report
    output_report = os.path.join(base_dir, "outputs", "ALS_GENOMIC_ATLAS.md")
    generate_report(db_path, output_report)
    
    # 6. Verification check for test run
    if test_run:
        logger.info("Running built-in gate check verification...")
        verify_db(db_path)
        logger.info("Verification PASSED!")

def verify_db(db_path: str):
    conn = duckdb.connect(db_path, read_only=True)
    try:
        # Check genes table
        ngenes = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        if ngenes == 0:
            raise ValueError("Verification failed: genes table is empty.")
        logger.info(f"Genes check passed: {ngenes} records.")

        # Check transcripts table
        ntranscripts = conn.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
        if ntranscripts == 0:
            raise ValueError("Verification failed: transcripts table is empty.")
        logger.info(f"Transcripts check passed: {ntranscripts} records.")

        # Check variants table
        nvariants = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        if nvariants == 0:
            raise ValueError("Verification failed: variants table is empty.")
        logger.info(f"Variants check passed: {nvariants} records.")

        # Check regulatory elements
        nreg = conn.execute("SELECT COUNT(*) FROM regulatory_elements").fetchone()[0]
        if nreg == 0:
            raise ValueError("Verification failed: regulatory_elements table is empty.")
        logger.info(f"Regulatory elements check passed: {nreg} records.")

        # Check expression
        nexpr = conn.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
        if nexpr == 0:
            raise ValueError("Verification failed: expression table is empty.")
        logger.info(f"Expression check passed: {nexpr} records.")

        # Check pathways_and_domains
        npw = conn.execute("SELECT COUNT(*) FROM pathways_and_domains").fetchone()[0]
        if npw == 0:
            raise ValueError("Verification failed: pathways_and_domains table is empty.")
        logger.info(f"Pathways & domains check passed: {npw} records.")

        # Check structures
        nstruct = conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
        if nstruct == 0:
            raise ValueError("Verification failed: structures table is empty.")
        logger.info(f"Structures check passed: {nstruct} records.")

        # Check foldseek_matches
        nfs = conn.execute("SELECT COUNT(*) FROM foldseek_matches").fetchone()[0]
        if nfs == 0:
            raise ValueError("Verification failed: foldseek_matches table is empty.")
        logger.info(f"Foldseek matches check passed: {nfs} records.")

        # Check foldseek_matched_drugs_trials
        nfs_drugs = conn.execute("SELECT COUNT(*) FROM foldseek_matched_drugs_trials").fetchone()[0]
        logger.info(f"Foldseek matched drugs check passed: {nfs_drugs} records.")

    finally:
        conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ALS Genomic Atlas Pipeline")
    parser.add_argument("--config", default=None, help="Path to configuration file")
    parser.add_argument("--test-run", action="store_true", help="Run a verification gate check with subset of genes")
    args = parser.parse_args()
    
    run_pipeline(args.config, args.test_run)
