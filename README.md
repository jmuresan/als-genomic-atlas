# ALS Genomic Atlas: Graph Data Dictionary

This document serves as a reference manual for the node and edge types, attributes, and biological data sources represented in the exported GraphML network (`outputs/als_genomic_atlas.graphml`).

<img width="1402" height="874" alt="Screenshot 2026-06-23 at 14 54 41" src="https://github.com/user-attachments/assets/01555743-6287-4322-a9cd-ecf46b79d0af" />
<img width="235" height="167" alt="Screenshot 2026-06-23 at 14 55 21" src="https://github.com/user-attachments/assets/e679076b-6033-470f-9c0e-eca42f6dc840" />


---

## Node Types

The network is a heterogeneous biological graph composed of **11 distinct node types** identified by the `node_type` attribute.

| Node Type (`node_type`) | Definition | Source Database(s) | Primary Attributes |
| :--- | :--- | :--- | :--- |
| **`gene`** | Represents human protein-coding genes. Includes 46 seed panel genes and their interacting partner genes. | Ensembl, UniProt, STRING | `ensembl_id` (string), `uniprot_id` (string), `chromosome` (string), `start_pos` (long), `end_pos` (long), `protein_description` (string), `als_seed` (long: `1` = seed panel, `0` = STRING network partner) |
| **`transcript`** | Represents mRNA transcript isoforms transcribed from a gene. | Ensembl | `mane_select` (boolean: True if MANE Select representative transcript), `length` (long), `exons` (long: exon count) |
| **`variant`** | Represents human genomic variants, primarily from ClinVar. | ClinVar, dbSNP, gnomAD, AlphaGenome | `clinical_significance` (string), `rsid` (string), `chromosome` (string), `position` (long), `hgvs` (string), `gnomad_pli` (double), `gnomad_loeuf` (double), `gnomad_allele_freq` (double), `alphagenome_consequence` (string), `alphagenome_pathogenicity` (double) |
| **`regulatory_element`** | Represents promoter or enhancer cis-regulatory elements associated with a gene. | ENCODE (SCREEN), UCSC Genome Browser | `element_type` (string: `'promoter'` or `'enhancer'`), `score` (double: ENCODE screening score), `ucsc_conservation_score` (double: phyloP/phastCons), `tfbs` (string: comma-separated list of transcription factor binding motifs) |
| **`tissue`** | Represents human tissues or anatomical locations. | GTEx Portal | None (label represents the tissue name, e.g., `'Brain - Spinal cord (cervical c-1)'`) |
| **`annotation`** | Represents biological pathways, Gene Ontology (GO) terms, or protein domain signatures. | Reactome, QuickGO (EBI), InterPro | `annotation_type` (string: `'pathway'`, `'go_term'`, or `'domain'`) |
| **`drug`** | Represents bioactive molecules, small compounds, or therapeutics. | Open Targets, ChEMBL | `max_clinical_phase` (double: e.g. `4.0` = Approved, `3.0` = Phase III), `mechanism_of_action` (string), `status` (string), `purpose` (string: disease indications) |
| **`trial`** | Represents clinical trials. | ClinicalTrials.gov, Open Targets | `status` (string: trial recruiting/completion status), `mechanism_of_action` (string), `max_clinical_phase` (double) |
| **`structure`** | Represents 3D protein structures, both predicted models and experimentally determined coordinates. | AlphaFold DB, RCSB PDB | `structure_type` (string: `'AlphaFold'` or `'PDB'`), `uniprot_id` (string), `plddt` (double: AlphaFold confidence score), `disorder_score` (double), `method` (string: PDB experimental method, e.g. `'X-RAY'`) |
| **`foldseek_target`** | Represents structural neighbor proteins (PDB or AlphaFold structures) showing high similarity. | Foldseek Search API | `db` (string: target structure database, e.g. `'afdb-swissprot'`) |
| **`disease`** | Represents diseases, syndromes, or clinical phenotypes. | ClinVar | `protein_description` (string: metadata for reference node) |

---

## Edge Types

Relationships between nodes are represented by **14 edge types** identified by the `edge_type` attribute. They are grouped into standard sub-networks using the `network` attribute.

| Edge Type (`edge_type`) | Source Node | Target Node | Sub-Network (`network`) | Definition / Source | Key Attributes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`has_transcript`** | `gene` | `transcript` | `gene-transcript` | Links a gene to its transcript isoforms. Source: **Ensembl** | None |
| **`has_variant`** | `gene` | `variant` | `gene-variant` | Links a variant to its host gene. Source: **ClinVar** | None |
| **`associated_with_disease`** | `variant` | `disease` | `variant-disease` | Links a variant to its ClinVar phenotype/disease classification. Source: **ClinVar** | None |
| **`associated_with_disease_direct`** | `gene` | `disease` | `gene-disease` | Aggregated direct relationship between a gene and a ClinVar phenotype based on variant counts. Source: **ClinVar** | `Weight` (double: variant count), `source_type` (string: `'clinvar'`) |
| **`seed_panel_disease`** | `gene` | `disease` | `gene-disease` | Connects genes in the seed panel that lack ClinVar variants directly to the reference phenotype. Source: **Seed Curation** | `Weight` (double: constant `1.0`), `source_type` (string: `'seed_panel'`) |
| **`regulated_by`** | `gene` | `regulatory_element` | `gene-regulatory` | Connects a gene to its regulating enhancer/promoter elements. Source: **ENCODE** | None |
| **`expressed_in`** | `gene` | `tissue` | `gene-bodyimpact` | Connects a gene to a tissue where its expression has been quantified. Source: **GTEx, HPA** | `tpm` (double: expression in TPM), `hpa_localization` (string: subcellular localization), `hpa_score` (string) |
| **`associated_with`** | `gene` | `annotation` | `gene-pathway` | Links a gene to a functional annotation (pathway, GO term, domain). Source: **Reactome, QuickGO, InterPro** | None |
| **`interacts_with`** | `gene` | `gene` | `string` | Represents physical or functional protein-protein interactions (sorted lexicographically). Source: **STRING DB** | `confidence_score` (double), `Weight` (double: identical to confidence score) |
| **`targeted_by`** | `gene` or `foldseek_target` | `drug` | `gene-drug-trial` or `foldseek-drug-trial` | Represents a compound targeting a seed gene or a Foldseek match target. Source: **Open Targets, ChEMBL** | `max_clinical_phase` (double), `mechanism_of_action` (string), `status` (string), `purpose` (string) |
| **`involved_in`** | `gene` or `foldseek_target` | `trial` | `gene-drug-trial` or `foldseek-drug-trial` | Represents a clinical trial targeting a seed gene or a Foldseek match target. Source: **ClinicalTrials.gov** | `max_clinical_phase` (double), `mechanism_of_action` (string), `status` (string), `purpose` (string) |
| **`has_structure`** | `gene` | `structure` | `gene-structure` | Connects a gene to its PDB or AlphaFold structures. Source: **AlphaFold DB, RCSB PDB** | None |
| **`structurally_similar_to`** | `gene` | `foldseek_target` | `foldseek` | Connects a query gene's AlphaFold structure model to structural neighbors. Source: **Foldseek** | `db` (string), `probability` (double), `query_coverage` (double), `evalue` (double), `seq_identity` (double), `alignment_length` (long) |
| **`similar_to`** | `drug` | `drug` | `drug-similarity` | Links an original drug to chemically similar repurposing candidate compounds. Source: **ChEMBL Similarity** | `similarity` (double: percent structural similarity), `max_clinical_phase` (double), `purpose` (string) |
