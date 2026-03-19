"""Curated pathway, disease, frequency, and literature evidence retrieval."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import requests

from prostr_alfo.config import Settings
from prostr_alfo.models.schemas import DiseaseAssociation, LiteratureRecord, Mutation, PathwayRecord, VariantEvidence

UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
REACTOME_PATHWAYS_URL = "https://reactome.org/ContentService/data/mapping/UniProt/{accession}/pathways"
PROTEINS_VARIATION_URL = "https://www.ebi.ac.uk/proteins/api/variation/{accession}"
PUBMED_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def _cache_path(settings: Settings, namespace: str, name: str) -> Path:
    target = settings.cache_dir / namespace
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def _fetch_json(url: str, cache_path: Path, params: dict[str, Any] | None = None) -> Any:
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    response = requests.get(url, params=params, headers={"Accept": "application/json"}, timeout=60)
    response.raise_for_status()
    data = response.json()
    cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def fetch_uniprot_record(accession: str, settings: Settings) -> dict[str, Any]:
    """Fetch UniProt metadata for one accession."""

    return _fetch_json(
        UNIPROT_JSON_URL.format(accession=accession),
        _cache_path(settings, "uniprot", f"{accession}.json"),
    )


def fetch_reactome_pathways(accession: str, settings: Settings) -> list[PathwayRecord]:
    """Fetch protein-level Reactome pathways."""

    data = _fetch_json(
        REACTOME_PATHWAYS_URL.format(accession=accession),
        _cache_path(settings, "reactome", f"{accession}.json"),
    )
    pathways: list[PathwayRecord] = []
    for item in data[:10]:
        stable_id = item.get("stId") or item.get("stIdVersion") or str(item.get("dbId"))
        pathways.append(
            PathwayRecord(
                source="Reactome",
                pathway_id=stable_id,
                name=item.get("displayName", "Unnamed pathway"),
                url=f"https://reactome.org/content/detail/{stable_id}",
            )
        )
    return pathways


def fetch_variation_dataset(accession: str, settings: Settings) -> dict[str, Any]:
    """Fetch curated protein variant records from the EBI Proteins API."""

    return _fetch_json(
        PROTEINS_VARIATION_URL.format(accession=accession),
        _cache_path(settings, "variation", f"{accession}.json"),
    )


def match_variant_feature(dataset: dict[str, Any], mutation: Mutation) -> dict[str, Any] | None:
    """Return the exact curated feature matching one amino-acid substitution."""

    for feature in dataset.get("features", []):
        if feature.get("type") != "VARIANT":
            continue
        if int(feature.get("begin", -1)) != mutation.position:
            continue
        if str(feature.get("wildType", "")).upper() != mutation.wild_type:
            continue
        mutated = str(feature.get("mutatedType") or feature.get("alternativeSequence") or "").upper()
        if mutated == mutation.mutant:
            return feature
    return None


def _clean_significances(feature: dict[str, Any]) -> list[str]:
    return sorted({entry.get("type", "").strip() for entry in feature.get("clinicalSignificances", []) if entry.get("type")})


def _frequency_summary(feature: dict[str, Any]) -> str | None:
    frequencies = []
    for entry in feature.get("populationFrequencies", []):
        name = entry.get("populationName") or entry.get("source") or "population"
        value = entry.get("frequency")
        if value is None:
            continue
        frequencies.append(f"{name}: {value}")
    if not frequencies:
        return None
    return "; ".join(frequencies[:5])


def _prediction_summary(feature: dict[str, Any]) -> list[str]:
    predictions = []
    for entry in feature.get("predictions", [])[:5]:
        algorithm = entry.get("predAlgorithmNameType", "Predictor")
        label = entry.get("predictionValType", "No label")
        score = entry.get("score")
        score_text = f" (score {score})" if score is not None else ""
        predictions.append(f"{algorithm}: {label}{score_text}")
    return predictions


def _disease_associations(feature: dict[str, Any]) -> list[DiseaseAssociation]:
    diseases: list[DiseaseAssociation] = []
    for association in feature.get("association", [])[:8]:
        db_references = [
            {
                "name": ref.get("name", ""),
                "id": ref.get("id", ""),
                "url": ref.get("url", ""),
            }
            for ref in association.get("dbReferences", [])
        ]
        evidence_pmids = [
            source.get("id", "")
            for evidence in association.get("evidences", [])
            for source in [evidence.get("source", {})]
            if source.get("name", "").lower() == "pubmed" and source.get("id")
        ]
        diseases.append(
            DiseaseAssociation(
                name=association.get("name", "Unnamed association"),
                description=association.get("description"),
                database_references=db_references,
                evidence_pmids=evidence_pmids,
            )
        )
    return diseases


def _source_links(feature: dict[str, Any]) -> list[dict[str, str]]:
    links = []
    for xref in feature.get("xrefs", [])[:10]:
        url = xref.get("url") or xref.get("alternativeUrl")
        if url:
            links.append({"name": xref.get("name", "Source"), "id": xref.get("id", ""), "url": url})
    return links


def fetch_pubmed_summaries(pmids: list[str], settings: Settings) -> list[LiteratureRecord]:
    """Fetch PubMed citation metadata for linked papers."""

    unique_pmids = sorted({pmid for pmid in pmids if pmid})[:10]
    if not unique_pmids:
        return []

    key = hashlib.sha256(",".join(unique_pmids).encode("utf-8")).hexdigest()[:16]
    data = _fetch_json(
        PUBMED_ESUMMARY_URL,
        _cache_path(settings, "pubmed", f"{key}.json"),
        params={"db": "pubmed", "id": ",".join(unique_pmids), "retmode": "json"},
    )
    result = data.get("result", {})
    citations: list[LiteratureRecord] = []
    for pmid in result.get("uids", []):
        entry = result.get(pmid, {})
        citations.append(
            LiteratureRecord(
                pmid=pmid,
                title=entry.get("title", "Untitled record"),
                journal=entry.get("fulljournalname"),
                publication_date=entry.get("pubdate"),
                authors=[author.get("name", "") for author in entry.get("authors", [])[:6] if author.get("name")],
                url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            )
        )
    return citations


def search_pubmed_variant_literature(
    *,
    gene_name: str | None,
    mutation: Mutation,
    settings: Settings,
) -> list[LiteratureRecord]:
    """Search PubMed for the exact variant string when curated papers are missing."""

    if not gene_name:
        return []

    query = f"{gene_name} {mutation.code} mutation"
    key = hashlib.sha256(query.encode("utf-8")).hexdigest()[:16]
    search_data = _fetch_json(
        PUBMED_ESEARCH_URL,
        _cache_path(settings, "pubmed", f"search_{key}.json"),
        params={"db": "pubmed", "term": query, "retmode": "json", "retmax": 3},
    )
    pmids = search_data.get("esearchresult", {}).get("idlist", [])
    return fetch_pubmed_summaries(pmids, settings)


def _phenotype_sufficiency(significances: list[str], matched: bool) -> tuple[str, str]:
    if not matched:
        return (
            "High uncertainty",
            "No exact curated record for this amino-acid substitution was found, so phenotype sufficiency is unknown.",
        )

    lower = {value.lower() for value in significances}
    if any("pathogenic" in value for value in lower) and not any("uncertain" in value for value in lower):
        return (
            "Moderate uncertainty",
            "Curated clinical databases support disease relevance for this variant, but phenotype sufficiency still depends on zygosity, penetrance, inheritance, and molecular context.",
        )
    if any("uncertain" in value for value in lower) or any("conflicting" in value for value in lower):
        return (
            "High uncertainty",
            "Curated sources classify this variant as uncertain or conflicting, so phenotype sufficiency is not established.",
        )
    return (
        "Moderate uncertainty",
        "No direct standalone phenotype claim was found; any disease consequence should be considered context-dependent.",
    )


def collect_variant_evidence(
    *,
    accession: str,
    mutations: list[Mutation],
    settings: Settings,
) -> list[VariantEvidence]:
    """Collect curated variant evidence for each mutation."""

    if not mutations:
        return []

    try:
        uniprot_record = fetch_uniprot_record(accession, settings)
        gene_name = (
            uniprot_record.get("genes", [{}])[0].get("geneName", {}).get("value")
            if uniprot_record.get("genes")
            else None
        )
        pathways = fetch_reactome_pathways(accession, settings)
        variation_dataset = fetch_variation_dataset(accession, settings)
    except Exception as exc:
        return [
            VariantEvidence(
                mutation=mutation,
                matched=False,
                uncertainty="High uncertainty",
                evidence_summary=f"External evidence lookup failed: {exc}",
                phenotype_sufficiency="No external clinical evidence could be retrieved.",
            )
            for mutation in mutations
        ]

    evidence_records: list[VariantEvidence] = []
    for mutation in mutations:
        feature = match_variant_feature(variation_dataset, mutation)
        if feature is None:
            literature = search_pubmed_variant_literature(gene_name=gene_name, mutation=mutation, settings=settings)
            uncertainty, phenotype_text = _phenotype_sufficiency([], matched=False)
            evidence_records.append(
                VariantEvidence(
                    mutation=mutation,
                    matched=False,
                    uncertainty=uncertainty,
                    evidence_summary=(
                        "No exact curated EBI Proteins variant record was found for this substitution. "
                        "Any linked papers come from an exact-string PubMed query and should be treated as supportive, not definitive."
                    ),
                    phenotype_sufficiency=phenotype_text,
                    pathways=pathways,
                    literature=literature,
                )
            )
            continue

        clinical_significances = _clean_significances(feature)
        diseases = _disease_associations(feature)
        curated_pmids = [pmid for disease in diseases for pmid in disease.evidence_pmids]
        literature = fetch_pubmed_summaries(curated_pmids, settings)
        if not literature:
            literature = search_pubmed_variant_literature(gene_name=gene_name, mutation=mutation, settings=settings)

        uncertainty, phenotype_text = _phenotype_sufficiency(clinical_significances, matched=True)
        evidence_records.append(
            VariantEvidence(
                mutation=mutation,
                matched=True,
                uncertainty=uncertainty,
                evidence_summary=(
                    "Exact curated variant evidence was found in the EBI Proteins variation API. "
                    "Linked papers come from curated disease-association records and can include gene- or phenotype-level evidence rather than a dedicated experiment for the exact substitution."
                ),
                clinical_significances=clinical_significances,
                consequence_type=feature.get("consequenceType"),
                frequency_summary=_frequency_summary(feature),
                phenotype_sufficiency=phenotype_text,
                pathways=pathways,
                diseases=diseases,
                literature=literature,
                source_links=_source_links(feature),
                prediction_summary=_prediction_summary(feature),
            )
        )
    return evidence_records
