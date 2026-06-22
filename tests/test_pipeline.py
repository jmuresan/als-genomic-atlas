import os
import tempfile
import yaml
import pytest
import duckdb
from src.config import Config
from src.ingest.cache import DiskCache
from src.db.schema import create_tables
from src.db.populate import populate_all
from src.atlas.generate_atlas import generate_report

@pytest.fixture
def temp_config_file():
    config_data = {
        "seed_genes": ["SOD1", "FUS"],
        "api_settings": {
            "offline_mode": True,
            "cache_dir": "data/raw/cache",
            "string_db": {
                "confidence_threshold": 0.7,
                "partner_limit": 10
            },
            "pubmed": {
                "limit_per_gene": 10
            }
        },
        "scoring_weights": {
            "open_targets_association": 0.25,
            "clinvar_pathogenicity": 0.20,
            "pathway_centrality": 0.15,
            "string_centrality": 0.15,
            "literature_volume": 0.15,
            "druggability": 0.10
        }
    }
    with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
        yaml.dump(config_data, f)
        temp_path = f.name
    yield temp_path
    os.remove(temp_path)

def test_config_parsing(temp_config_file):
    cfg = Config(temp_config_file)
    assert cfg.seed_genes == ["SOD1", "FUS"]
    assert cfg.offline_mode is True
    assert cfg.string_partner_limit == 10
    
    # Boundary validation check
    with pytest.raises(ValueError):
        cfg.data["api_settings"]["string_db"]["partner_limit"] = -1
        # re-trigger init logic with bad data
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as bad_f:
            yaml.dump(cfg.data, bad_f)
            bad_path = bad_f.name
        try:
            Config(bad_path)
        finally:
            os.remove(bad_path)

def test_cache_mechanism():
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = DiskCache(tmpdir, offline_mode=False)
        key = cache.generate_cache_key("test_source", "test_endpoint", {"param": "value"})
        assert len(key) == 64
        
        # Test read/write
        test_data = {"result": "success"}
        cache.write("test_source", "test_endpoint", {"param": "value"}, test_data)
        
        cached_val = cache.read("test_source", "test_endpoint", {"param": "value"})
        assert cached_val == test_data

def test_database_population():
    # Setup in-memory DuckDB
    conn = duckdb.connect(":memory:")
    create_tables(conn)
    
    # Assert tables are created
    tables = [t[0] for t in conn.execute("SHOW TABLES").fetchall()]
    assert "genes" in tables
    assert "variants" in tables
    assert "transcripts" in tables
    assert "regulatory_elements" in tables
    assert "expression" in tables
    assert "pathways_and_domains" in tables
    assert "structures" in tables
    
    # Mock data package
    mock_data = {
        "gene": "SOD1",
        "uniprot_id": "P00441",
        "ensembl": {
            "ensembl_id": "ENSG00000142168",
            "transcripts": [
                {"id": "ENST00000270142", "mane_select": True, "length": 2000, "exons": 5}
            ],
            "coordinates": {"chr": "21", "start": 31659693, "end": 31668931}
        },
        "uniprot": {
            "primaryAccession": "P00441",
            "genes": [{"geneName": {"value": "SOD1"}}],
            "proteinDescription": {"recommendedName": {"fullName": {"value": "Superoxide dismutase"}}}
        },
        "clinvar": {
            "result": {
                "uids": ["12345"],
                "12345": {
                    "uid": "12345",
                    "accession": "VCV000012345",
                    "germline_classification": {
                        "description": "Pathogenic",
                        "trait_set": [{"trait_name": "Amyotrophic lateral sclerosis"}]
                    }
                }
            }
        },
        "dbsnp": {"rsid": "rs121912442", "chromosome": "21", "position": 31668406, "hgvs": "NC_C>T"},
        "gnomad": {"pli": 0.12, "loeuf": 0.85, "allele_freq": 0.0001},
        "alphagenome": {"non_coding_variant_consequence": "regulatory_disruption", "pathogenicity_score": 0.8},
        "encode": {
            "promoters": [{"id": "EH38E1", "score": 0.9}],
            "enhancers": []
        },
        "ucsc": {"conservation_score": 0.9, "tfbs": ["MA0139.1"]},
        "gtex": {"tissues": [{"tissue": "Spinal cord", "tpm": 120.0}]},
        "hpa": {"localization": "Cytoplasm", "score": "High"},
        "reactome": [{"stId": "R-HSA-70326", "displayName": "Superoxide degradation"}],
        "quickgo": {"results": [{"goId": "GO:0004784", "goName": "superoxide dismutase activity"}]},
        "interpro": {"results": [{"metadata": {"accession": "IPR001424", "name": "Superoxide dismutase, copper/zinc binding site"}}]},
        "string": [{"preferredName_A": "SOD1", "preferredName_B": "CCS", "score": 0.99}],
        "open_targets": {
            "approvedSymbol": "SOD1",
            "drugs": [{"maxClinicalStage": "PHASE_III", "drug": {"id": "CHEMBL1201484", "name": "Riluzole", "maximumClinicalStage": "APPROVAL"}}]
        },
        "clinical_trials": {"trials": [{"nct_id": "NCT00000123", "title": "Riluzole Trial", "status": "COMPLETED"}]},
        "alphafold": {"plddt": 95.0, "disorder_score": 0.05},
        "pdb": {"pdb_ids": ["1HLN"], "method": "X-RAY"}
    }
    
    populate_all(conn, mock_data)
    
    # Assert database populated correctly
    ngenes = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    assert ngenes == 1
    
    ntx = conn.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
    assert ntx == 1
    
    nvar = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    assert nvar == 1
    
    conn.close()

def test_report_generation():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test.duckdb")
        report_path = os.path.join(tmpdir, "test_report.md")
        
        conn = duckdb.connect(db_path)
        create_tables(conn)
        
        # Insert minimal seed gene for report
        conn.execute("""
            INSERT INTO genes (gene_symbol, ensembl_id, uniprot_id, chromosome, start_pos, end_pos, protein_description)
            VALUES ('SOD1', 'ENSG00000142168', 'P00441', '21', 31659693, 31668931, 'Superoxide dismutase')
        """)
        conn.close()
        
        generate_report(db_path, report_path)
        
        assert os.path.exists(report_path)
        with open(report_path, "r", encoding="utf-8") as f:
            content = f.read()
            assert "# ALS Genomic Atlas" in content
            assert "## SOD1" in content
