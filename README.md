# ALS Genomic Atlas

Batch pipeline that queries public bioinformatics APIs for a fixed panel of amyotrophic lateral sclerosis (ALS)–associated human genes, normalizes responses into DuckDB, and renders a single Markdown reference document.

Primary artifact: `outputs/ALS_GENOMIC_ATLAS.md` (one section per gene, eleven annotation categories).

This repository does not redistribute raw third-party databases. It stores API responses in a local cache and derived tables in DuckDB. Redistribution or publication of extracted data requires compliance with each provider’s license and citation requirements (see Credits).

## Gene panel

Forty-six gene symbols in `config.yaml` (`seed_genes`):

- Forty-one human genes annotated with UniProt keyword KW-0036 (“Amyotrophic lateral sclerosis”).
- Five additional genes from the OMIM ALS phenotypic series (PS105400) and GeneReviews ALS literature: CFAP410, MOBP, SCFD1, TAF15, UNC13A.

OMIM and GeneReviews informed gene selection only; this pipeline does not query OMIM as a data source.

## Data categories and API sources

| # | Category | Sources (as implemented in `src/ingest/client.py`) |
|---|----------|-----------------------------------------------------|
| 1 | Gene and transcript mapping | Ensembl REST, NCBI Entrez (nuccore), UniProt REST |
| 2 | Variants and pathogenicity | ClinVar and dbSNP (NCBI E-utilities), gnomAD GraphQL API, AlphaGenome API |
| 3 | Transcriptional regulation and epigenomics | ENCODE SCREEN / Factorbook GraphQL, UCSC hub API, JASPAR REST, UniBind REST |
| 4 | Expression and tissue specificity | GTEx Portal API v2, Human Protein Atlas download API |
| 5 | Pathways and functional annotation | Reactome Content Service, EBI QuickGO, InterPro REST |
| 6 | Protein–protein interactions | STRING REST |
| 7 | Clinical translation and druggability | Open Targets Platform GraphQL, ChEMBL REST, ClinicalTrials.gov API v2 |
| 8 | Three-dimensional structures | AlphaFold Database API, RCSB PDB search API |
| 9 | Structural similarity | Foldseek web API (`search.foldseek.com`) |
| 10 | Matched-target drugs and trials | Open Targets (targets resolved from Foldseek hits) |
| 11 | Similar compounds and repurposing | ChEMBL similarity search, Open Targets |

Categories 9–11 are sequential: Foldseek searches use AlphaFold-derived coordinates; downstream drug and trial fields depend on Open Targets and ChEMBL.

## Pipeline behavior

`run_pipeline.py` executes in order:

1. Load and validate `config.yaml`.
2. Create or update DuckDB schema (`src/db/schema.py`).
3. For each seed gene, `IngestionManager` fetches all categories and `populate_all` writes normalized rows (`src/db/populate.py`).
4. `generate_report` renders Markdown from DuckDB (`src/atlas/generate_atlas.py`).

HTTP responses are cached under `data/raw/cache/` keyed by SHA-256 of the request (`src/ingest/cache.py`). With `api_settings.offline_mode: true`, cache misses raise an error instead of calling the network.

On request failure, individual clients return empty structures; the run continues for other genes and endpoints. AlphaGenome, UCSC conservation, and UniBind endpoints are often unavailable and may yield empty sections.

NCBI E-utilities calls in this codebase do not currently set the recommended `tool` and `email` query parameters. Operators should add them before high-volume or production use (see NCBI policy in Credits).

## Requirements

- Python 3.12
- Dependencies: `requirements.txt`

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

Full panel:

```bash
python run_pipeline.py
```

First five genes with post-run table checks:

```bash
python run_pipeline.py --test-run
```

Alternate config path:

```bash
python run_pipeline.py --config path/to/config.yaml
```

## Outputs

| Path | Description |
|------|-------------|
| `data/processed/als_genomic_atlas.duckdb` | Normalized relational tables |
| `outputs/ALS_GENOMIC_ATLAS.md` | Rendered atlas |
| `data/raw/cache/` | Per-request JSON cache (regenerable; not versioned in git) |

## Configuration

`config.yaml` fields used by the pipeline:

- `seed_genes` — gene symbols to process.
- `api_settings.offline_mode` — network vs cache-only.
- `api_settings.cache_dir` — cache directory.
- `api_settings.string_db` — STRING confidence threshold and partner limit.
- `api_settings.foldseek` — probability threshold and maximum hits.
- `scoring_weights` — weights for internal priority scoring (normalized to sum 1.0 at load).

## Tests

```bash
pytest
```

## Repository layout

```
run_pipeline.py
config.yaml
src/config.py
src/ingest/client.py
src/ingest/cache.py
src/db/schema.py
src/db/populate.py
src/atlas/generate_atlas.py
tests/test_pipeline.py
video/als_atlas_eli5.py   # optional Manim explainer; not part of the data pipeline
```

## Credits and attribution

The atlas aggregates facts retrieved from the resources below. If you publish work that uses this pipeline or its outputs, cite the underlying databases and tools according to each provider’s policy. Wording below follows public “cite us” / terms pages as of June 2026; confirm current requirements on each site before publication.

### Ensembl

- **Site:** https://www.ensembl.org  
- **Attribution:** Cite the Ensembl overview article for the release you used and note the Ensembl release version. Current overview: Dyer SC, et al. Ensembl 2025. *Nucleic Acids Res.* 53(D1):D948–D957.  
- **Policy:** https://www.ensembl.org/info/about/legal/

### NCBI / National Library of Medicine

- **Sites:** https://www.ncbi.nlm.nih.gov (Entrez, ClinVar, dbSNP, nuccore)  
- **Attribution:** Acknowledge the National Library of Medicine (NLM) and National Center for Biotechnology Information (NCBI) as the source. NLM-authored material is in the public domain; downstream use should include appropriate NLM acknowledgment.  
- **E-utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25497/ — identify your application with `tool` and `email` parameters; observe usage frequency limits.  
- **Molecular databases:** https://www.ncbi.nlm.nih.gov/home/about/policies/

### ClinVar

- **Site:** https://www.ncbi.nlm.nih.gov/clinvar/  
- **Attribution:** Attribute ClinVar as a data source; cite a ClinVar publication (e.g. Landrum MJ, et al. ClinVar: improvements to accessing data. *Nucleic Acids Res.* — PMID 29165669).  
- **Use:** https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/

### UniProt

- **Site:** https://www.uniprot.org  
- **Attribution:** The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* 53:D609–D617 (2025). For API use, cite Ahmad S, et al. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* (2025) (REST API).  
- **License:** Creative Commons Attribution 4.0 (CC BY 4.0) for copyrightable database content — https://rest.uniprot.org/help/license

### gnomAD

- **Site:** https://gnomad.broadinstitute.org  
- **Attribution:** Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581, 434–443 (2020). PMID 32461654.  
- **License:** Data distributed under CC BY 4.0 (confirm on https://gnomad.broadinstitute.org/faq).

### AlphaGenome

- **Site:** https://www.alphagenome.ai / https://api.alphagenome.org  
- **Attribution:** Follow citation and terms published by the AlphaGenome provider for the API version and model release you query.

### ENCODE, SCREEN, and Factorbook

- **Sites:** https://www.encodeproject.org , https://screen.wenglab.org , Factorbook GraphQL (`factorbook.api.wenglab.org`)  
- **Attribution:** Cite ENCODE Consortium publications, acknowledge the producing laboratory, and reference ENCODE accessions (ENCSR…, ENCFF…) where applicable. See https://www.encodeproject.org/help/citing-encode/

### UCSC Genome Browser

- **Site:** https://genome.ucsc.edu  
- **Attribution:** Cite the Genome Browser update article (see https://genome.ucsc.edu/cite.html) and reference http://genome.ucsc.edu.  
- **Graphics / data:** https://genome.ucsc.edu/license/

### JASPAR

- **Site:** https://jaspar.elixir.no  
- **Attribution:** Khan A, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. *Nucleic Acids Res.* (2024).  
- **Policy:** https://jaspar.elixir.no/cite

### UniBind

- **Site:** https://unibind.uio.no  
- **Attribution:** Pinter RC, et al. UniBind–a database of direct TF-DNA interactions. *Nucleic Acids Res.* 51(D1):D185–D191 (2023).  
- **License:** CC BY 4.0 (stated on provider site).

### GTEx Portal

- **Site:** https://gtexportal.org  
- **Attribution:** GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369, 1318–1330 (2020). PMID 32913098.

### Human Protein Atlas

- **Site:** https://www.proteinatlas.org  
- **Attribution:** Uhlén M, et al. Proteomics. Tissue-based map of the human proteome. *Science* 347, 1260419 (2015). PMID 28392948.  
- **License:** CC BY-SA 4.0 — https://www.proteinatlas.org/about/licence

### Reactome

- **Site:** https://reactome.org  
- **Attribution:** Gillespie M, et al. Reactome 2022. *Nucleic Acids Res.* (2022) and related pathway citations listed at https://reactome.org/cite

### Gene Ontology and QuickGO

- **Sites:** https://geneontology.org , https://www.ebi.ac.uk/QuickGO/  
- **Attribution:** Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. *Nucleic Acids Res.* (2023); Ashburner M, et al. Gene ontology: tool for the unification of biology. *Nat Genet.* 25, 25–29 (2000). QuickGO: Binns D, et al. QuickGO: a web-based tool for Gene Ontology searching. *Bioinformatics* (2009).  
- **License:** CC BY 4.0 — https://geneontology.org/docs/go-citation-policy/

### InterPro

- **Site:** https://www.ebi.ac.uk/interpro/  
- **Attribution:** Cite the InterPro database article listed at https://www.ebi.ac.uk/interpro/ and EBI terms of use: https://www.ebi.ac.uk/about/terms-of-use/

### STRING

- **Site:** https://string-db.org  
- **Attribution:** Szklarczyk D, et al. STRING v12.0: updated protein–protein association networks with enriched annotation coverage across diverse organisms. *Nucleic Acids Res.* (2023).  
- **References:** https://string-db.org/cgi/about?footer_active_subpage=references

### Open Targets Platform

- **Site:** https://platform.opentargets.org  
- **Attribution:** Buniello A, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Res.* (2025).  
- **License:** Platform data CC0 1.0; code Apache 2.0 — https://platform-docs.opentargets.org/licence.md

### ChEMBL

- **Site:** https://www.ebi.ac.uk/chembl/  
- **Attribution:** Bento AP, et al. ChEMBL web services: streamlining access to drug discovery data and utilities. *Nucleic Acids Res.* (2023).  
- **License:** CC BY-SA 3.0 (ChEMBL data; confirm in ChEMBL FAQ).

### ClinicalTrials.gov

- **Site:** https://clinicaltrials.gov  
- **Attribution:** Acknowledge ClinicalTrials.gov and the National Library of Medicine as the source of trial metadata.  
- **Terms:** https://clinicaltrials.gov/about-site/terms-conditions

### AlphaFold Protein Structure Database

- **Site:** https://alphafold.ebi.ac.uk  
- **Attribution:** Bertoni D, et al. AlphaFold Protein Structure Database 2025. *Nucleic Acids Res.* (2025). For model methodology: Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021).  
- **License:** CC BY 4.0 for structure files — https://alphafold.ebi.ac.uk/download

### RCSB Protein Data Bank

- **Site:** https://www.rcsb.org  
- **Attribution:** Cite the PDB entry (PDB ID), the primary structure publication, RCSB PDB, and the molecular viewer used (e.g. Mol*).  
- **Policy:** https://www.rcsb.org/pages/policies

### Foldseek

- **Site:** https://search.foldseek.com (API); https://github.com/steineggerlab/foldseek  
- **Attribution:** van Kempen M, et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol.* (2024).

---

Gene panel curation references (selection criteria, not ingested as APIs): UniProt keyword KW-0036; OMIM phenotypic series PS105400; GeneReviews ALS.