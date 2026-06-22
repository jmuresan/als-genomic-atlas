# ALS Genomic Atlas

A data pipeline that pulls multi-omics annotations for a panel of ALS-associated
genes from public bioinformatics APIs, normalizes them into a DuckDB database, and
renders a single Markdown reference document.

The output is [`outputs/ALS_GENOMIC_ATLAS.md`](outputs/ALS_GENOMIC_ATLAS.md): one
section per gene, organized into 11 categories from transcript structure through
variant pathogenicity, regulation, expression, pathways, interactions, drugs, and
3D structure.

## Gene panel

46 genes (`seed_genes` in [`config.yaml`](config.yaml)): the 41 human genes carrying
the UniProt "Amyotrophic lateral sclerosis" keyword (KW-0036), plus CFAP410, MOBP,
SCFD1, TAF15, and UNC13A from the OMIM ALS phenotypic series (PS105400) and
GeneReviews.

## Data categories and sources

| # | Category | Sources |
|---|----------|---------|
| 1 | Gene & transcript mapping | Ensembl, NCBI nuccore, UniProt |
| 2 | Variants & pathogenicity | ClinVar, dbSNP, gnomAD, AlphaGenome |
| 3 | Transcriptional regulation & epigenomics | ENCODE SCREEN/factorbook, UCSC conservation, JASPAR, UniBind |
| 4 | Expression & tissue specificity | GTEx, Human Protein Atlas |
| 5 | Pathways & functional annotation | Reactome, QuickGO (Gene Ontology), InterPro |
| 6 | Network interactions | STRING |
| 7 | Clinical translation & druggability | Open Targets, ChEMBL, ClinicalTrials.gov |
| 8 | 3D structural biology | AlphaFold, RCSB PDB |
| 9 | Structural similarity | Foldseek (queried with the AlphaFold model) |
| 10 | Matched-target drugs & trials | Foldseek hits resolved through Open Targets |
| 11 | Similar compounds & repurposing | ChEMBL similarity, Open Targets |

Categories 9–11 chain off each other: Foldseek finds structural homologs of each
gene's predicted structure, those targets are looked up in Open Targets for known
drugs and trials, and ChEMBL similarity search expands those drugs into repurposing
candidates.

## How it works

`run_pipeline.py` runs the stages in order:

1. Load [`config.yaml`](config.yaml) (gene list, API settings, scoring weights).
2. Create the DuckDB schema ([`src/db/schema.py`](src/db/schema.py)).
3. For each gene, fetch all categories through `IngestionManager`
   ([`src/ingest/client.py`](src/ingest/client.py)) and write the normalized rows
   into DuckDB ([`src/db/populate.py`](src/db/populate.py)).
4. Render the Markdown atlas from the database
   ([`src/atlas/generate_atlas.py`](src/atlas/generate_atlas.py)).

Every HTTP response is cached to disk under `data/raw/cache/` as a JSON file keyed
by a SHA-256 hash of the request ([`src/ingest/cache.py`](src/ingest/cache.py)).
Reruns serve from cache instead of hitting the network. A failed API call returns
an empty record rather than aborting the run, so one flaky endpoint does not sink
the whole panel.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Requires Python 3.12.

## Usage

Run the full panel:

```bash
python run_pipeline.py
```

Quick gate check on the first 5 genes, with table-population assertions at the end:

```bash
python run_pipeline.py --test-run
```

Point at a different config:

```bash
python run_pipeline.py --config path/to/config.yaml
```

### Offline mode

Set `api_settings.offline_mode: true` in the config to serve everything from the
cache. A cache miss in this mode raises an error instead of falling back to the
network, which is useful for reproducible reruns.

## Outputs

| Path | Contents |
|------|----------|
| `data/processed/als_genomic_atlas.duckdb` | Normalized tables (genes, variants, structures, foldseek matches, …) |
| `outputs/ALS_GENOMIC_ATLAS.md` | The rendered atlas, one section per gene |
| `data/raw/cache/` | Raw API responses, one JSON file per request |

The cache is regenerable and is not tracked in git (see
[`.gitignore`](.gitignore)). It rebuilds from the public APIs on the next run.

## Configuration

Key fields in [`config.yaml`](config.yaml):

- `seed_genes` — the gene panel to process.
- `api_settings.offline_mode` — serve from cache only.
- `api_settings.cache_dir` — where cached responses live.
- `api_settings.string_db` — STRING confidence threshold and partner limit.
- `api_settings.foldseek` — probability threshold and hit cap for structural matches.
- `scoring_weights` — relative weights for the gene-priority score; normalized to
  sum to 1.0 at load.

## Tests

```bash
pytest
```

Covers config parsing and the populate/render path against cached data.

## Layout

```
run_pipeline.py            Entry point and orchestration
config.yaml                Gene panel and API/scoring settings
src/config.py              Config loader and validation
src/ingest/client.py       API clients and the IngestionManager
src/ingest/cache.py        SHA-256 keyed disk cache
src/db/schema.py           DuckDB table definitions
src/db/populate.py         Maps API responses into tables
src/atlas/generate_atlas.py  Renders the Markdown atlas
tests/test_pipeline.py     Tests
```

## Notes

Data comes from public sources with their own licenses and citation requirements;
check each source before redistributing. AlphaGenome, UCSC conservation, and
UniBind clients depend on endpoints that are frequently unavailable and will return
empty records when they fail rather than stopping the run.
