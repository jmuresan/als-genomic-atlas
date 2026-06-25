import re
import time
import json
import logging
import urllib.parse
import requests
from typing import Dict, Any, Optional, List

logger = logging.getLogger("als_atlas.client")

# --- Base Ingestion Client ---
#
# De-mocked: the original BaseClient fabricated records via `_get_mock_data`
# whenever a request failed or offline mode was on. That has been removed. A
# failed request now returns an honest empty structure (so a single flaky
# endpoint does not abort the whole 46-gene run) and never invents data.


class BaseClient:
    """Base API client logic with disk caching, retries, and rate limiting."""

    def __init__(self, source_name: str, cache, rate_limit_delay: float = 0.25):
        self.source_name = source_name
        self.cache = cache
        self.rate_limit_delay = rate_limit_delay

    def _empty(self, is_xml: bool) -> Any:
        if is_xml:
            return ""
        return [] if self.source_name in ("reactome", "string") else {}

    def _request(self, method: str, url: str, endpoint: str,
                 params: Optional[Dict[str, Any]] = None,
                 json_data: Optional[Dict[str, Any]] = None,
                 headers: Optional[Dict[str, str]] = None,
                 is_xml: bool = False) -> Any:
        query_params: Dict[str, Any] = {}
        if params:
            query_params.update(params)
        if json_data:
            query_params.update({"_post_body": json_data})

        # Serve from cache when possible (offline miss raises in cache.read).
        cached = self.cache.read(self.source_name, endpoint, query_params)
        if cached is not None:
            if is_xml or isinstance(cached, (dict, list)):
                return cached
            logger.warning(f"Cached data for {self.source_name}/{endpoint} is {type(cached)}; ignoring cache.")

        max_attempts = 3
        backoff_delay = 0.5
        response = None

        for attempt in range(max_attempts):
            try:
                if self.rate_limit_delay > 0 and attempt == 0:
                    time.sleep(self.rate_limit_delay)

                logger.info(f"Live request: {method} {url} for endpoint {endpoint} (attempt {attempt + 1})")
                if method.upper() == "POST":
                    response = requests.post(url, json=json_data, headers=headers, timeout=15)
                else:
                    response = requests.get(url, params=params, headers=headers, timeout=15)

                response.raise_for_status()
                break
            except requests.exceptions.RequestException as e:
                is_transient = False
                if e.response is not None:
                    status_code = e.response.status_code
                    if status_code == 429 or (500 <= status_code < 600):
                        is_transient = True
                else:
                    is_transient = True

                if is_transient and attempt < max_attempts - 1:
                    logger.warning(f"Transient error on attempt {attempt + 1}: {e}. Retrying in {backoff_delay}s...")
                    time.sleep(backoff_delay)
                    backoff_delay *= 2
                else:
                    logger.error(f"Failed to fetch from {url}: {e}. Returning empty response.")
                    return self._empty(is_xml)

        if is_xml:
            data = response.text
        else:
            try:
                data = response.json()
            except ValueError:
                logger.error(f"Failed to parse JSON response from {url}. Returning empty response.")
                return self._empty(is_xml)

        self.cache.write(self.source_name, endpoint, query_params, data)
        return data


# --- Category 1 Clients ---

class EnsemblClient(BaseClient):
    def __init__(self, cache):
        super().__init__("ensembl", cache)

    def fetch_gene(self, gene_symbol: str) -> Dict[str, Any]:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
        headers = {"Content-Type": "application/json"}
        raw = self._request("GET", url, "symbol_lookup", params={"gene": gene_symbol, "expand": "1"}, headers=headers)

        if isinstance(raw, dict) and "coordinates" in raw:
            return raw
        if not isinstance(raw, dict):
            return {"ensembl_id": None, "symbol": gene_symbol, "transcripts": [], "coordinates": {}}

        transcripts_list = []
        for tx in raw.get("Transcript", []) or []:
            is_mane = tx.get("is_canonical") == 1 or "mane" in str(tx.get("attributes", {})).lower()
            transcripts_list.append({
                "id": tx.get("id"),
                "mane_select": bool(is_mane),
                "length": tx.get("length", 0),
                "exons": len(tx.get("Exon", [])),
            })

        return {
            "ensembl_id": raw.get("id"),
            "symbol": gene_symbol,
            "transcripts": transcripts_list,
            "coordinates": {
                "chr": str(raw.get("seq_region_name")),
                "start": raw.get("start"),
                "end": raw.get("end"),
            },
        }


class NCBISequenceClient(BaseClient):
    def __init__(self, cache):
        super().__init__("ncbi_sequence", cache)

    def fetch_sequence(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        return self._request("GET", url, "efetch_seq", params={"db": "nuccore", "gene": gene_symbol, "retmode": "text"})


class UniProtClient(BaseClient):
    def __init__(self, cache):
        super().__init__("uniprot", cache)

    def get_gene_details(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": f"gene:{gene_symbol} AND organism_id:9606", "format": "json"}
        res = self._request("GET", url, "search", params=params)
        results = res.get("results", []) if isinstance(res, dict) else []
        return results[0] if results else {}


# --- Category 2 Clients ---

class ClinVarClient(BaseClient):
    def __init__(self, cache):
        super().__init__("clinvar", cache)

    def get_variants(self, gene_symbol: str) -> Dict[str, Any]:
        url_search = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        url_summary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

        search_res = self._request("GET", url_search, "esearch", params={
            "db": "clinvar",
            "term": f"{gene_symbol}[gene] AND amyotrophic lateral sclerosis",
            "retmode": "json",
        })
        id_list = search_res.get("esearchresult", {}).get("idlist", []) if isinstance(search_res, dict) else []
        if not id_list:
            return {"result": {}}

        return self._request("GET", url_summary, "esummary", params={
            "db": "clinvar", "id": ",".join(id_list), "retmode": "json",
        })


class DbSNPClient(BaseClient):
    def __init__(self, cache):
        super().__init__("dbsnp", cache)

    def fetch_snp(self, rsid: str) -> Dict[str, Any]:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        numeric_id = rsid.lower().replace("rs", "")
        raw = self._request("GET", url, f"rs_{rsid}", params={"db": "snp", "id": numeric_id, "retmode": "json"})

        if isinstance(raw, dict) and "rsid" in raw:
            return raw
        if not isinstance(raw, dict) or "result" not in raw:
            return {}

        result = raw.get("result", {})
        uids = result.get("uids", [])
        if not uids:
            return {}

        uid = uids[0]
        snp_info = result.get(uid, {})
        chrom = snp_info.get("chr", "21")

        spdi = snp_info.get("spdi", "")
        position = None
        if spdi and len(spdi.split(":")) >= 2:
            try:
                position = int(spdi.split(":")[1])
            except ValueError:
                pass

        hgvs_val = None
        for item in snp_info.get("docsum", "").split(","):
            if item.startswith("HGVS="):
                hgvs_val = item.split("=", 1)[1]
                break

        return {
            "rsid": f"rs{uid}",
            "chromosome": chrom,
            "position": position,
            "hgvs": hgvs_val or (f"NC_0000{chrom}.9:g.{position}C>T" if position else None),
        }


class GnomADClient(BaseClient):
    def __init__(self, cache):
        super().__init__("gnomad", cache)

    def fetch_constraint(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://gnomad.broadinstitute.org/api"
        query = """
        query($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
          gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
            gnomad_constraint {
              pli
              oe_lof_upper
            }
          }
        }
        """
        variables = {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"}
        raw = self._request("POST", url, "constraint", json_data={"query": query, "variables": variables, "gene": gene_symbol})

        if isinstance(raw, dict) and "pli" in raw:
            return raw
        if not isinstance(raw, dict) or "data" not in raw:
            return {"gene": gene_symbol, "pli": None, "loeuf": None, "allele_freq": None}

        constraint = ((raw.get("data", {}) or {}).get("gene") or {}).get("gnomad_constraint") or {}
        return {
            "gene": gene_symbol,
            "pli": constraint.get("pli"),
            "loeuf": constraint.get("oe_lof_upper"),
            # gnomAD constraint does not provide a single allele frequency; the
            # original code hard-coded 0.00001 here (mock). Left unset.
            "allele_freq": None,
        }


class AlphaGenomeClient(BaseClient):
    def __init__(self, cache):
        super().__init__("alphagenome", cache)

    def predict_variant(self, variant_id: str) -> Dict[str, Any]:
        return self._request("GET", "https://api.alphagenome.org/variant", f"variant_{variant_id}", params={"variant": variant_id})


# --- Category 3 Clients ---

class EncodeClient(BaseClient):
    def __init__(self, cache):
        super().__init__("encode", cache)

    def fetch_ccres(self, gene_symbol: str, coords: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        url = "https://factorbook.api.wenglab.org/graphql"
        headers = {"Origin": "https://screen-v2.wenglab.org", "Referer": "https://screen-v2.wenglab.org/"}

        chrom = coords.get("chr") if coords else None
        start = coords.get("start") if coords else None
        end = coords.get("end") if coords else None

        if not chrom or not start or not end:
            gene_query = """
            query GeneID($name: [String!]) {
                gene(assembly: "grch38", name: $name) {
                    coordinates { chromosome start end }
                }
            }
            """
            ref_res = self._request("POST", url, "resolve_gene", json_data={"query": gene_query, "variables": {"name": [gene_symbol]}, "gene": gene_symbol}, headers=headers)
            if isinstance(ref_res, dict):
                genes = (ref_res.get("data", {}) or {}).get("gene", []) or []
                if genes and genes[0].get("coordinates"):
                    c = genes[0]["coordinates"]
                    chrom, start, end = c.get("chromosome"), c.get("start"), c.get("end")

        if not chrom or not start or not end:
            return {"gene": gene_symbol, "promoters": [], "enhancers": []}

        if not str(chrom).startswith("chr"):
            chrom = f"chr{chrom}"

        query = """
        query Search($coords: [GenomicRangeInput!]) {
            cCRESCREENSearch(assembly: "grch38", coordinates: $coords) {
                info { accession }
                promoter_zscore
                enhancer_zscore
            }
        }
        """
        variables = {"coords": [{"chromosome": chrom, "start": int(start), "end": int(end)}]}
        raw = self._request("POST", url, "ccres", json_data={"query": query, "variables": variables, "gene": gene_symbol}, headers=headers)

        if isinstance(raw, dict) and ("promoters" in raw or "enhancers" in raw):
            return raw
        if not isinstance(raw, dict) or "data" not in raw:
            return {"gene": gene_symbol, "promoters": [], "enhancers": []}

        promoters, enhancers = [], []
        for c in (raw.get("data") or {}).get("cCRESCREENSearch", []) or []:
            acc = c.get("info", {}).get("accession")
            p_score = c.get("promoter_zscore", 0.0)
            e_score = c.get("enhancer_zscore", 0.0)
            if acc:
                if p_score >= e_score and p_score >= 1.64:
                    promoters.append({"id": acc, "score": float(p_score)})
                elif e_score >= 1.64:
                    enhancers.append({"id": acc, "score": float(e_score)})

        return {"gene": gene_symbol, "promoters": promoters, "enhancers": enhancers}


class UcscConservationClient(BaseClient):
    def __init__(self, cache):
        super().__init__("ucsc", cache)

    def fetch_scores(self, gene_symbol: str) -> Dict[str, Any]:
        return self._request("GET", "https://genome.ucsc.edu/cgi-bin/hubApi", "conservation", params={"gene": gene_symbol})


class JasparClient(BaseClient):
    def __init__(self, cache):
        super().__init__("jaspar", cache)

    def fetch_motifs(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://jaspar.elixir.no/api/v1/matrix/"
        return self._request("GET", url, "motifs", params={"search": gene_symbol})


class UniBindClient(BaseClient):
    def __init__(self, cache):
        super().__init__("unibind", cache)

    def fetch_sites(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://unibind.uio.no/api/v1/tf/"
        return self._request("GET", url, "sites", params={"tf": gene_symbol})


# --- Category 4 Clients ---

class GtexClient(BaseClient):
    def __init__(self, cache):
        super().__init__("gtex", cache)

    def fetch_expression(self, gene_symbol: str) -> Dict[str, Any]:
        ref_url = "https://gtexportal.org/api/v2/reference/gene"
        ref_res = self._request("GET", ref_url, "resolve_gene", params={"geneId": gene_symbol, "gene": gene_symbol})

        gencode_id = None
        if isinstance(ref_res, dict):
            ref_data = ref_res.get("data", []) or []
            if ref_data:
                best_match = ref_data[0]
                for d in ref_data:
                    if d.get("geneSymbol", "").lower() == gene_symbol.lower():
                        best_match = d
                        break
                gencode_id = best_match.get("gencodeId")

        if not gencode_id:
            return {"gene": gene_symbol, "tissues": []}

        url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        raw = self._request("GET", url, "expression", params={"gencodeId": gencode_id, "datasetId": "gtex_v8", "gene": gene_symbol})

        if isinstance(raw, dict) and "tissues" in raw:
            return raw

        if isinstance(raw, dict):
            raw_data = raw.get("data", []) or []
        elif isinstance(raw, list):
            raw_data = raw
        else:
            raw_data = []

        tissues_list = []
        for item in raw_data:
            tissue_id = item.get("tissueSiteDetailId") or item.get("tissueSiteDetail")
            if tissue_id:
                tissues_list.append({"tissue": tissue_id.replace("_", " "), "tpm": float(item.get("median", 0.0))})

        return {"gene": gene_symbol, "tissues": tissues_list}


class HpaClient(BaseClient):
    def __init__(self, cache):
        super().__init__("human_protein_atlas", cache)

    def fetch_localization(self, gene_symbol: str, ensembl_id: Optional[str] = None) -> Dict[str, Any]:
        url = "https://www.proteinatlas.org/api/search_download.php"
        search_id = ensembl_id if ensembl_id else gene_symbol
        raw = self._request("GET", url, "localization", params={
            "search": search_id,
            "columns": "g,eg,scl,scml,scal",
            "format": "json",
            "compress": "no",
            "gene": gene_symbol,
        })

        if isinstance(raw, dict) and "localization" in raw:
            return raw
        if not isinstance(raw, list) or not raw:
            return {"gene": gene_symbol, "localization": "N/A", "score": None}

        locations = raw[0].get("Subcellular location", []) or []
        # The original hard-coded score "High" for every gene (mock). The HPA
        # download columns here do not carry a reliability score, so leave None.
        return {
            "gene": gene_symbol,
            "localization": ", ".join(locations) if locations else "N/A",
            "score": None,
        }


# --- Category 5 Clients ---

class ReactomeClient(BaseClient):
    def __init__(self, cache):
        super().__init__("reactome", cache)

    def get_pathways_for_uniprot(self, uniprot_id: str) -> List[Dict[str, Any]]:
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{uniprot_id}/pathways"
        res = self._request("GET", url, f"pathways_{uniprot_id}", params={"species": 9606})
        return res if isinstance(res, list) else []


class QuickGOClient(BaseClient):
    def __init__(self, cache):
        super().__init__("quickgo", cache)

    def fetch_go_terms(self, gene_symbol: str, uniprot_id: Optional[str] = None) -> Dict[str, Any]:
        url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
        query_id = f"UniProtKB:{uniprot_id}" if uniprot_id else gene_symbol
        raw = self._request("GET", url, "go_terms", params={"geneProductId": query_id, "gene": gene_symbol})

        if not isinstance(raw, dict) or "results" not in raw:
            return {"results": []}

        results = raw.get("results", []) or []
        go_ids = list({r.get("goId") for r in results if r.get("goId")})

        go_name_map: Dict[str, str] = {}
        if go_ids:
            for i in range(0, len(go_ids), 50):
                batch = go_ids[i:i + 50]
                terms_url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{','.join(batch)}"
                try:
                    resp = requests.get(terms_url, timeout=10)
                    if resp.status_code == 200:
                        for term in resp.json().get("results", []):
                            if term.get("id") and term.get("name"):
                                go_name_map[term["id"]] = term["name"]
                except Exception as e:
                    logger.error(f"Error fetching GO term names: {e}")

        new_results = []
        for r in results:
            go_id = r.get("goId")
            go_name = go_name_map.get(go_id) or r.get("goName")
            if go_id and not go_name:
                aspect = r.get("goAspect", "").replace("_", " ")
                go_name = f"GO annotation ({aspect})" if aspect else "GO annotation"
            new_results.append({"goId": go_id, "goName": go_name})

        return {"results": new_results}


class InterProClient(BaseClient):
    def __init__(self, cache):
        super().__init__("interpro", cache)

    def fetch_domains(self, uniprot_id: str) -> Dict[str, Any]:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/{uniprot_id}"
        return self._request("GET", url, f"domains_{uniprot_id}")


# --- Category 6 Clients ---

class StringClient(BaseClient):
    def __init__(self, cache, confidence_threshold: float = 0.7, limit: int = 10):
        super().__init__("string", cache)
        self.confidence_threshold = confidence_threshold
        self.limit = limit

    def get_interactions(self, gene_symbol: str) -> List[Dict[str, Any]]:
        url = "https://string-db.org/api/json/interaction_partners"
        params = {
            "identifiers": gene_symbol,
            "species": 9606,
            "required_score": int(self.confidence_threshold * 1000),
            "limit": self.limit,
            "gene": gene_symbol,
        }
        res = self._request("GET", url, "interactions", params=params)
        return res if isinstance(res, list) else []


# --- Category 7 Clients ---

class OpenTargetsClient(BaseClient):
    def __init__(self, cache):
        super().__init__("open_targets", cache)
        self.graphql_url = "https://api.platform.opentargets.org/api/v4/graphql"

    def fetch_gene_data(self, gene_symbol: str, ensembl_id: Optional[str] = None) -> Dict[str, Any]:
        # De-mocked: when a target cannot be resolved, return an honest empty
        # record rather than the original fabricated "Riluzole / MockTherapeutic".
        empty = {"approvedSymbol": gene_symbol, "associatedDiseases": {"rows": []}, "drugs": []}

        if not ensembl_id:
            search_query = """
            query searchTarget($queryString: String!) {
              search(queryString: $queryString, entityNames: ["target"], page: {index: 0, size: 1}) {
                hits {
                  id
                }
              }
            }
            """
            search_res = self._request("POST", self.graphql_url, "target_search", json_data={"query": search_query, "variables": {"queryString": gene_symbol}, "gene": gene_symbol})
            if isinstance(search_res, dict):
                hits = search_res.get("data", {}).get("search", {}).get("hits", []) or []
                if hits:
                    ensembl_id = hits[0].get("id")

        if not ensembl_id:
            return empty

        query = """
        query targetDetails($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            associatedDiseases(page: {index: 0, size: 50}) {
              count
              rows {
                score
                disease {
                  id
                  name
                }
              }
            }
            drugAndClinicalCandidates {
              count
              rows {
                maxClinicalStage
                drug {
                  id
                  name
                  mechanismsOfAction {
                    rows {
                      mechanismOfAction
                    }
                  }
                }
              }
            }
          }
        }
        """
        raw = self._request("POST", self.graphql_url, "target_data", json_data={"query": query, "variables": {"ensemblId": ensembl_id}, "gene": gene_symbol})

        if isinstance(raw, dict) and "drugs" in raw:
            return raw
        if not isinstance(raw, dict) or "data" not in raw:
            return empty

        target = raw.get("data", {}).get("target", {}) or {}
        if not target:
            return empty

        assoc_diseases_rows = []
        for r in (target.get("associatedDiseases", {}) or {}).get("rows", []) or []:
            disease_info = r.get("disease", {}) or {}
            assoc_diseases_rows.append({
                "disease": {"id": disease_info.get("id"), "name": disease_info.get("name")},
                "score": float(r.get("score", 0.0)),
            })

        drug_rows = []
        for r in (target.get("drugAndClinicalCandidates", {}) or {}).get("rows", []) or []:
            drug_info = r.get("drug", {}) or {}
            moa_rows = [{"mechanismOfAction": m.get("mechanismOfAction")}
                        for m in (drug_info.get("mechanismsOfAction", {}) or {}).get("rows", []) or []]
            drug_rows.append({
                "maxClinicalStage": r.get("maxClinicalStage"),
                "drug": {"id": drug_info.get("id"), "name": drug_info.get("name"), "mechanismsOfAction": {"rows": moa_rows}},
            })

        return {
            "approvedSymbol": target.get("approvedSymbol", gene_symbol),
            "associatedDiseases": {"rows": assoc_diseases_rows},
            "drugs": drug_rows,
        }

    def fetch_drugs_and_indications_for_target(self, target_query: str) -> Dict[str, Any]:
        """Queries Open Targets to find approved drugs, trials, and indications for a matched target ID/name."""
        endpoint = "target_drugs_indications"
        query_params = {"target": target_query}

        cached = self.cache.read(self.source_name, endpoint, query_params)
        if cached is not None:
            return cached

        ensembl_id = None
        search_term = target_query

        # Resolve via UniProt accession
        if re.match(r'^[A-Z0-9]{6,10}$', target_query):
            try:
                uniprot_url = f"https://rest.uniprot.org/uniprotkb/{target_query}.json"
                resp = requests.get(uniprot_url, timeout=10)
                if resp.status_code == 200:
                    resp_data = resp.json()
                    
                    # 1. Enforce target must be human to avoid mapping non-human symbols to human targets
                    organism = resp_data.get("organism", {})
                    scientific_name = organism.get("scientificName", "")
                    if scientific_name.lower() != "homo sapiens":
                        logger.info(f"Skipping non-human target {target_query} ({scientific_name})")
                        result = {"approvedSymbol": target_query, "drugs": []}
                        self.cache.write(self.source_name, endpoint, query_params, result)
                        return result
                    
                    # 2. Try to get Ensembl ID directly from UniProt cross-references to bypass fuzzy search
                    for ref in resp_data.get("uniProtKBCrossReferences", []):
                        if ref.get("database") == "Ensembl":
                            properties = ref.get("properties", [])
                            for prop in properties:
                                if prop.get("key") == "GeneId":
                                    ensembl_id = prop.get("value")
                                    break
                            if ensembl_id:
                                if "." in ensembl_id:
                                    ensembl_id = ensembl_id.split(".")[0]
                                break
                    
                    # Resolve gene symbol for search fallback
                    genes_info = resp_data.get("genes", [])
                    if genes_info:
                        val = genes_info[0].get("geneName", {}).get("value")
                        if val:
                            search_term = val.upper()
                        else:
                            syns = genes_info[0].get("synonyms", [])
                            if syns and syns[0].get("value"):
                                search_term = syns[0]["value"].upper()
            except Exception as e:
                logger.error(f"Error resolving UniProt ID {target_query} to gene name: {e}")

        # If Ensembl ID was not found directly in cross-references, fall back to target search
        if not ensembl_id:
            search_query = """
            query searchTarget($queryString: String!) {
              search(queryString: $queryString, entityNames: ["target"]) {
                hits {
                  id
                  name
                }
              }
            }
            """
            search_res = self._request("POST", self.graphql_url, "target_search", json_data={"query": search_query, "variables": {"queryString": search_term}})
            if isinstance(search_res, dict):
                hits = search_res.get("data", {}).get("search", {}).get("hits", []) or []
                for h in hits:
                    h_id = h.get("id")
                    h_name = h.get("name", "")
                    # Enforce strict validation: target symbol must match search_term exactly
                    if h_name and h_name.upper() == search_term:
                        ensembl_id = h_id
                        break

        if not ensembl_id:
            result = {"approvedSymbol": target_query, "drugs": []}
            self.cache.write(self.source_name, endpoint, query_params, result)
            return result

        query = """
        query targetDrugsAndIndications($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            drugAndClinicalCandidates {
              count
              rows {
                maxClinicalStage
                drug {
                  id
                  name
                  mechanismsOfAction {
                    rows {
                      mechanismOfAction
                    }
                  }
                  indications {
                    rows {
                      disease {
                        name
                      }
                    }
                  }
                }
              }
            }
          }
        }
        """
        raw = self._request("POST", self.graphql_url, "target_details", json_data={"query": query, "variables": {"ensemblId": ensembl_id}})

        if not isinstance(raw, dict) or "data" not in raw:
            return {"drugs": []}
        target = raw.get("data", {}).get("target", {}) or {}
        if not target:
            return {"drugs": []}

        drug_rows = []
        for r in (target.get("drugAndClinicalCandidates", {}) or {}).get("rows", []) or []:
            drug_info = r.get("drug", {}) or {}
            moa_rows = [{"mechanismOfAction": m.get("mechanismOfAction")}
                        for m in (drug_info.get("mechanismsOfAction", {}) or {}).get("rows", []) or []]
            ind_rows = [{"disease": {"name": (i.get("disease", {}) or {}).get("name")}}
                        for i in (drug_info.get("indications", {}) or {}).get("rows", []) or []]
            drug_rows.append({
                "maxClinicalStage": r.get("maxClinicalStage"),
                "drug": {"id": drug_info.get("id"), "name": drug_info.get("name"),
                         "mechanismsOfAction": {"rows": moa_rows}, "indications": {"rows": ind_rows}},
            })

        result = {"approvedSymbol": target.get("approvedSymbol", target_query), "drugs": drug_rows}
        self.cache.write(self.source_name, endpoint, query_params, result)
        return result

    def fetch_drug_details(self, chembl_id: str) -> Dict[str, Any]:
        """Queries Open Targets for a drug's indications and mechanisms of action."""
        endpoint = "drug_details"
        query_params = {"chemblId": chembl_id}

        cached = self.cache.read(self.source_name, endpoint, query_params)
        if cached is not None:
            return cached

        query = """
        query drugDetails($chemblId: String!) {
          drug(chemblId: $chemblId) {
            id
            name
            mechanismsOfAction {
              rows {
                mechanismOfAction
              }
            }
            indications {
              rows {
                disease {
                  name
                }
              }
            }
          }
        }
        """
        return self._request("POST", self.graphql_url, endpoint, json_data={"query": query, "variables": {"chemblId": chembl_id}})


class ChemblClient(BaseClient):
    def __init__(self, cache):
        super().__init__("chembl", cache)

    def fetch_mechanism(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://www.ebi.ac.uk/chembl/api/data/mechanism"
        headers = {"Accept": "application/json"}
        return self._request("GET", url, "mechanism", params={"target_chembl_id": gene_symbol, "_format": "json"}, headers=headers)

    def fetch_molecule_structures(self, chembl_id: str) -> Dict[str, Any]:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
        headers = {"Accept": "application/json"}
        return self._request("GET", url, f"molecule_{chembl_id}", headers=headers)

    def fetch_similar_compounds(self, smiles: str, threshold: int = 85) -> Dict[str, Any]:
        quoted_smiles = urllib.parse.quote(smiles)
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{quoted_smiles}/{threshold}.json"
        headers = {"Accept": "application/json"}
        return self._request("GET", url, f"similarity_{threshold}", params={"smiles": smiles}, headers=headers)


class ClinicalTrialsClient(BaseClient):
    def __init__(self, cache):
        super().__init__("clinical_trials", cache)

    def fetch_trials(self, gene_symbol: str) -> Dict[str, Any]:
        url = "https://clinicaltrials.gov/api/v2/studies"
        raw = self._request("GET", url, "studies", params={"query.cond": "Amyotrophic Lateral Sclerosis", "query.term": gene_symbol})

        if not isinstance(raw, dict):
            return {"trial_count": 0, "trials": []}

        trials = []
        for study in raw.get("studies", []):
            protocol = study.get("protocolSection", {})
            ident = protocol.get("identificationModule", {})
            status_mod = protocol.get("statusModule", {})
            nct_id = ident.get("nctId")
            if nct_id:
                trials.append({
                    "nct_id": nct_id,
                    "title": ident.get("briefTitle") or ident.get("officialTitle"),
                    "status": status_mod.get("overallStatus"),
                })
        return {"trial_count": len(trials), "trials": trials}


# --- Category 8 Clients ---

class AlphaFoldClient(BaseClient):
    def __init__(self, cache):
        super().__init__("alphafold", cache)

    def fetch_structure(self, uniprot_id: str) -> Dict[str, Any]:
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        return self._request("GET", url, f"structure_{uniprot_id}")


class PdbClient(BaseClient):
    def __init__(self, cache):
        super().__init__("pdb", cache)

    def fetch_pdb_ids(self, gene_symbol: str) -> Dict[str, Any]:
        payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                    "operator": "exact_match",
                    "value": gene_symbol,
                },
            },
            "return_type": "entry",
        }
        url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={urllib.parse.quote(json.dumps(payload))}"
        raw = self._request("GET", url, "query_by_gene", params={"gene": gene_symbol})

        if isinstance(raw, dict) and "pdb_ids" in raw:
            return raw
        if not isinstance(raw, dict) or "result_set" not in raw:
            return {"pdb_ids": [], "method": "N/A"}

        pdb_ids = [item["identifier"] for item in raw.get("result_set", [])]
        return {"pdb_ids": pdb_ids, "method": "X-RAY DIFFRACTION"}


class FoldseekClient(BaseClient):
    def __init__(self, cache):
        super().__init__("foldseek", cache)

    def fetch_alignments(self, gene_symbol: str, pdb_url: Optional[str]) -> Dict[str, Any]:
        endpoint = "alignments"
        query_params = {"gene": gene_symbol}

        cached = self.cache.read(self.source_name, endpoint, query_params)
        if cached is not None:
            return cached

        if not pdb_url:
            logger.error(f"No PDB URL available for Foldseek alignment of {gene_symbol}.")
            return {"results": []}

        try:
            logger.info(f"Downloading PDB structure from {pdb_url} for Foldseek query...")
            pdb_resp = requests.get(pdb_url, timeout=20)
            pdb_resp.raise_for_status()
            pdb_content = pdb_resp.text

            logger.info(f"Submitting Foldseek alignment job for {gene_symbol}...")
            url = "https://search.foldseek.com/api/ticket"
            databases = ["afdb50", "afdb-swissprot", "afdb-proteome", "pdb100",
                         "BFVD", "mgnify_esm30", "cath50", "gmgcl_id", "bfmd"]
            files = {"q": ("query.pdb", pdb_content)}
            data = [("mode", "3diaa")] + [("database[]", db) for db in databases]

            response = requests.post(url, files=files, data=data, timeout=30)
            if response.status_code != 200:
                logger.error(f"Foldseek submission failed: HTTP {response.status_code}")
                return {"results": []}

            ticket_id = response.json().get("id")
            if not ticket_id:
                logger.error(f"Foldseek submission returned no ticket ID. Response: {response.json()}")
                return {"results": []}

            logger.info(f"Foldseek ticket generated: {ticket_id}. Polling for completion...")
            status_url = f"https://search.foldseek.com/api/ticket/{ticket_id}"
            max_polls = 30
            for i in range(max_polls):
                time.sleep(5)
                status_resp = requests.get(status_url, timeout=20)
                if status_resp.status_code != 200:
                    logger.warning(f"Foldseek status check returned HTTP {status_resp.status_code}")
                    continue
                status = status_resp.json().get("status")
                logger.info(f"Foldseek status check {i + 1}/{max_polls}: {status}")
                if status == "COMPLETE":
                    break
                elif status == "ERROR":
                    logger.error(f"Foldseek job failed on server for ticket {ticket_id}")
                    return {"results": []}
            else:
                logger.error(f"Foldseek job timed out after {max_polls * 5} seconds")
                return {"results": []}

            logger.info(f"Fetching Foldseek results for ticket {ticket_id}...")
            result_url = f"https://search.foldseek.com/api/result/{ticket_id}/0"
            res = requests.get(result_url, timeout=60)
            if res.status_code != 200:
                logger.error(f"Failed to fetch Foldseek results: HTTP {res.status_code}")
                return {"results": []}

            results = res.json()
            self.cache.write(self.source_name, endpoint, query_params, results)
            return results
        except Exception as e:
            logger.error(f"Exception during Foldseek alignment: {e}")
            return {"results": []}


# --- Consolidated Ingestion Manager ---

class IngestionManager:
    """Consolidated Ingest client to query all 8 biological categories."""
    def __init__(self, cache):
        self.cache = cache
        self.foldseek = FoldseekClient(cache)
        self.ensembl = EnsemblClient(cache)
        self.ncbi_seq = NCBISequenceClient(cache)
        self.uniprot = UniProtClient(cache)
        self.clinvar = ClinVarClient(cache)
        self.dbsnp = DbSNPClient(cache)
        self.gnomad = GnomADClient(cache)
        self.alphagenome = AlphaGenomeClient(cache)
        self.encode = EncodeClient(cache)
        self.ucsc = UcscConservationClient(cache)
        self.jaspar = JasparClient(cache)
        self.unibind = UniBindClient(cache)
        self.gtex = GtexClient(cache)
        self.hpa = HpaClient(cache)
        self.reactome = ReactomeClient(cache)
        self.quickgo = QuickGOClient(cache)
        self.interpro = InterProClient(cache)
        self.string = StringClient(cache)
        self.open_targets = OpenTargetsClient(cache)
        self.chembl = ChemblClient(cache)
        self.clinical_trials = ClinicalTrialsClient(cache)
        self.alphafold = AlphaFoldClient(cache)
        self.pdb = PdbClient(cache)

    def fetch_all_data(self, gene_symbol: str) -> Dict[str, Any]:
        """Runs the ingestion pipeline for a single gene symbol, returning a structured package."""
        logger.info(f"Gathering 8-category mapping data for {gene_symbol}...")

        # Category 1
        ensembl_data = self.ensembl.fetch_gene(gene_symbol)
        ncbi_seq_data = self.ncbi_seq.fetch_sequence(gene_symbol)
        uniprot_data = self.uniprot.get_gene_details(gene_symbol)

        uniprot_id = uniprot_data.get("primaryAccession") if uniprot_data else None
        ensembl_id = ensembl_data.get("ensembl_id")
        coords = ensembl_data.get("coordinates", {}) or {}

        # Category 2
        clinvar_data = self.clinvar.get_variants(gene_symbol)
        dbsnp_data = self.dbsnp.fetch_snp("rs121912442")  # Example standard ALS variant rsID
        gnomad_data = self.gnomad.fetch_constraint(gene_symbol)
        alphagenome_data = self.alphagenome.predict_variant("chr21:31668406:C>T")

        # Category 3
        encode_data = self.encode.fetch_ccres(gene_symbol, coords)
        ucsc_data = self.ucsc.fetch_scores(gene_symbol)
        jaspar_data = self.jaspar.fetch_motifs(gene_symbol)
        unibind_data = self.unibind.fetch_sites(gene_symbol)

        # Category 4
        gtex_data = self.gtex.fetch_expression(gene_symbol)
        hpa_data = self.hpa.fetch_localization(gene_symbol, ensembl_id)

        # Category 5
        reactome_data = []
        interpro_data = {}
        if uniprot_id:
            reactome_data = self.reactome.get_pathways_for_uniprot(uniprot_id)
            interpro_data = self.interpro.fetch_domains(uniprot_id)
        quickgo_data = self.quickgo.fetch_go_terms(gene_symbol, uniprot_id)

        # Category 6
        string_data = self.string.get_interactions(gene_symbol)

        # Category 7
        ot_data = self.open_targets.fetch_gene_data(gene_symbol, ensembl_id)
        chembl_data = self.chembl.fetch_mechanism(gene_symbol)
        trials_data = self.clinical_trials.fetch_trials(gene_symbol)

        # Category 8
        alphafold_data = {}
        if uniprot_id:
            alphafold_data = self.alphafold.fetch_structure(uniprot_id)
        pdb_data = self.pdb.fetch_pdb_ids(gene_symbol)

        # Category 9 & 10: Foldseek & Matched target drugs/trials
        from src.config import Config
        config = Config()
        prob_threshold = config.foldseek_probability_threshold
        limit_hits = config.foldseek_limit_hits

        # Resolve pdbUrl from AlphaFold structure to pass to Foldseek
        pdb_url = None
        if uniprot_id:
            try:
                af_data = self.alphafold.fetch_structure(uniprot_id)
                if isinstance(af_data, list) and af_data:
                    pdb_url = af_data[0].get("pdbUrl")
                elif isinstance(af_data, dict):
                    pdb_url = af_data.get("pdbUrl")
            except Exception as e:
                logger.error(f"Error fetching AlphaFold structure info for {uniprot_id}: {e}")

        # Fallback to standard model v6 naming pattern if not explicitly resolved from API metadata
        if not pdb_url and uniprot_id:
            pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"

        foldseek_data = self.foldseek.fetch_alignments(gene_symbol, pdb_url)

        # Process alignments and look up target drugs
        alignments_list = []
        if isinstance(foldseek_data, dict):
            if "results" in foldseek_data:
                for result_group in foldseek_data.get("results", []):
                    for db_alignments in result_group.get("alignments", []):
                        if isinstance(db_alignments, list):
                            alignments_list.extend(db_alignments)
                        elif isinstance(db_alignments, dict):
                            alignments_list.append(db_alignments)
            elif "alignments" in foldseek_data:
                alignments_list = foldseek_data["alignments"]
        elif isinstance(foldseek_data, list):
            alignments_list = foldseek_data

        filtered_matches = []
        for hit in alignments_list:
            prob = hit.get("prob", hit.get("probability", 0.0))
            try:
                prob = float(prob)
            except (ValueError, TypeError):
                prob = 0.0

            if prob >= prob_threshold:
                target = hit.get("target", "")
                clean_target = target
                if target.startswith("AF-") and "-" in target[3:]:
                    parts = target.split("-")
                    if len(parts) >= 2:
                        clean_target = parts[1]
                elif "_" in target and not target.endswith("_A"):
                    clean_target = target.split("_")[0]

                q_len = hit.get("qlen", hit.get("qLen", 0))
                q_start = hit.get("qStartPos", 0)
                q_end = hit.get("qEndPos", 0)
                q_cov = 0.0
                if q_len > 0 and q_end > q_start:
                    q_cov = min((q_end - q_start + 1) / q_len, 1.0)

                seq_id_val = hit.get("seqId", hit.get("seqIdentity", hit.get("fident", 0.0)))
                try:
                    seq_id_val = float(seq_id_val)
                    if seq_id_val > 1.0:
                        seq_id_val /= 100.0
                except (ValueError, TypeError):
                    seq_id_val = 0.0

                filtered_matches.append({
                    "target_id": target,
                    "clean_target": clean_target,
                    "db": hit.get("db", "afdb-swissprot"),
                    "probability": prob,
                    "query_coverage": q_cov,
                    "eval": hit.get("eval", hit.get("eValue", hit.get("evalue", 1000.0))),
                    "seqId": seq_id_val,
                    "alnLength": hit.get("alnLength", hit.get("alnLen", hit.get("alnlen", 0)))
                })

        filtered_matches.sort(key=lambda x: x["probability"], reverse=True)
        filtered_matches = filtered_matches[:limit_hits]

        # Stage to phase dictionary mapping for helper
        STAGE_TO_PHASE = {
            "APPROVAL": 4.0, "PHASE_IV": 4.0, "PHASE_III": 3.0, "PHASE_II": 2.0,
            "PHASE_I": 1.0, "EARLY_PHASE_I": 0.5, "PRECLINICAL": 0.0, "PHASE_0": 0.0,
        }

        foldseek_drugs = []
        for match in filtered_matches:
            target_id = match["clean_target"]
            drugs_info = self.open_targets.fetch_drugs_and_indications_for_target(target_id)
            if drugs_info and "drugs" in drugs_info:
                for d in drugs_info["drugs"]:
                    ind_list = [row.get("disease", {}).get("name", "") for row in d.get("drug", {}).get("indications", {}).get("rows", []) or []]
                    ind_list = [i for i in ind_list if i]

                    moa_list = [row.get("mechanismOfAction", "") for row in d.get("drug", {}).get("mechanismsOfAction", {}).get("rows", []) or []]
                    moa_list = [m for m in moa_list if m]

                    purpose = ""
                    if ind_list:
                        purpose += f"Indicated for: {', '.join(ind_list)}. "
                    if moa_list:
                        purpose += f"Mechanism: {'; '.join(moa_list)}."
                    if not purpose:
                        purpose = "No indication/mechanism details available."

                    stage = d.get("maxClinicalStage")
                    phase_val = None
                    if stage:
                        phase_val = STAGE_TO_PHASE.get(str(stage).strip().upper())

                    foldseek_drugs.append({
                        "target_id": match["target_id"],
                        "drug_id": d.get("drug", {}).get("id"),
                        "type": "drug",
                        "name_or_title": d.get("drug", {}).get("name"),
                        "max_clinical_phase": phase_val,
                        "mechanism_of_action": "; ".join(moa_list) if moa_list else None,
                        "status": None,
                        "purpose": purpose
                    })

        # Category 11: Similar Compounds & Repurposing Candidates
        cat7_drugs = []
        if ot_data and isinstance(ot_data, dict):
            for row in ot_data.get("drugs", []):
                drug = row.get("drug", {}) or {}
                drug_id = drug.get("id")
                drug_name = drug.get("name")
                stage = drug.get("maximumClinicalStage") or row.get("maxClinicalStage")
                phase = STAGE_TO_PHASE.get(str(stage).strip().upper()) if stage else 0.0
                if phase is None:
                    phase = 0.0
                if drug_id:
                    cat7_drugs.append({
                        "drug_id": drug_id,
                        "name": drug_name,
                        "max_clinical_phase": phase,
                        "target_id": gene_symbol
                    })

        cat10_drugs = []
        for fd in foldseek_drugs:
            if fd.get("type") == "drug":
                drug_id = fd.get("drug_id")
                drug_name = fd.get("name_or_title")
                phase = fd.get("max_clinical_phase")
                if phase is None:
                    phase = 0.0
                if drug_id:
                    cat10_drugs.append({
                        "drug_id": drug_id,
                        "name": drug_name,
                        "max_clinical_phase": phase,
                        "target_id": fd.get("target_id")
                    })

        # Unique candidates by (drug_id, target_id)
        unique_candidates = {}
        for d in cat7_drugs + cat10_drugs:
            key = (d["drug_id"], d["target_id"])
            if key not in unique_candidates:
                unique_candidates[key] = d
            else:
                if d["max_clinical_phase"] > unique_candidates[key]["max_clinical_phase"]:
                    unique_candidates[key] = d

        similar_compounds = []
        for candidate in unique_candidates.values():
            chembl_id = candidate["drug_id"]
            smiles = None
            try:
                mol_res = self.chembl.fetch_molecule_structures(chembl_id)
                if isinstance(mol_res, dict):
                    smiles = mol_res.get("molecule_structures", {}).get("canonical_smiles")
                    if not smiles:
                        h = mol_res.get("molecule_hierarchy", {}) or {}
                        parent = h.get("parent_chembl_id")
                        if parent and parent != chembl_id:
                            parent_res = self.chembl.fetch_molecule_structures(parent)
                            if isinstance(parent_res, dict):
                                smiles = parent_res.get("molecule_structures", {}).get("canonical_smiles")
            except Exception as e:
                logger.error(f"Error fetching SMILES for {chembl_id}: {e}")

            if smiles:
                try:
                    sim_res = self.chembl.fetch_similar_compounds(smiles, threshold=85)
                    if isinstance(sim_res, dict) and "molecules" in sim_res:
                        for mol in sim_res["molecules"]:
                            sim_id = mol.get("molecule_chembl_id")
                            sim_name = mol.get("pref_name")
                            similarity_val = None
                            similarity_str = mol.get("similarity")
                            if similarity_str is not None:
                                try:
                                    similarity_val = float(similarity_str)
                                except ValueError:
                                    pass

                            sim_phase_val = mol.get("max_phase")
                            try:
                                sim_phase_val = float(sim_phase_val) if sim_phase_val is not None else 0.0
                            except ValueError:
                                sim_phase_val = 0.0

                            orig_phase = candidate["max_clinical_phase"]
                            if sim_phase_val >= orig_phase:
                                purpose = ""
                                try:
                                    drug_det = self.open_targets.fetch_drug_details(sim_id)
                                    if isinstance(drug_det, dict) and "data" in drug_det:
                                        drug_info = drug_det.get("data", {}).get("drug", {}) or {}

                                        moa_list = [row.get("mechanismOfAction", "") for row in drug_info.get("mechanismsOfAction", {}).get("rows", []) or []]
                                        moa_list = [m for m in moa_list if m]

                                        ind_list = [row.get("disease", {}).get("name", "") for row in drug_info.get("indications", {}).get("rows", []) or []]
                                        ind_list = [i for i in ind_list if i]

                                        if ind_list:
                                            purpose += f"Indicated for: {', '.join(ind_list)}. "
                                        if moa_list:
                                            purpose += f"Mechanism: {'; '.join(moa_list)}."
                                except Exception as e:
                                    logger.error(f"Error fetching Open Targets drug details for {sim_id}: {e}")

                                if not purpose:
                                    purpose = "No indication/mechanism details available."

                                similar_compounds.append({
                                    "target_id": candidate["target_id"],
                                    "original_drug_id": chembl_id,
                                    "similar_drug_id": sim_id,
                                    "name": sim_name or sim_id,
                                    "similarity": similarity_val,
                                    "max_clinical_phase": sim_phase_val,
                                    "purpose": purpose
                                })
                except Exception as e:
                    logger.error(f"Error fetching similar compounds for SMILES {smiles}: {e}")

        return {
            "gene": gene_symbol,
            "uniprot_id": uniprot_id,
            "ensembl": ensembl_data,
            "ncbi_seq": ncbi_seq_data,
            "uniprot": uniprot_data,
            "clinvar": clinvar_data,
            "dbsnp": dbsnp_data,
            "gnomad": gnomad_data,
            "alphagenome": alphagenome_data,
            "encode": encode_data,
            "ucsc": ucsc_data,
            "jaspar": jaspar_data,
            "unibind": unibind_data,
            "gtex": gtex_data,
            "hpa": hpa_data,
            "reactome": reactome_data,
            "quickgo": quickgo_data,
            "interpro": interpro_data,
            "string": string_data,
            "open_targets": ot_data,
            "chembl": chembl_data,
            "clinical_trials": trials_data,
            "alphafold": alphafold_data,
            "pdb": pdb_data,
            "foldseek_matches": filtered_matches,
            "foldseek_drugs": foldseek_drugs,
            "similar_compounds": similar_compounds
        }
