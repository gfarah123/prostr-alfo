"""Shared end-to-end analysis pipeline for CLI and Streamlit."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from prostr_alfo.analysis.evidence import collect_variant_evidence
from prostr_alfo.analysis.features import analyze_structure
from prostr_alfo.analysis.mutation import assess_mutations
from prostr_alfo.config import Settings, get_settings
from prostr_alfo.io.input import resolve_protein_input
from prostr_alfo.models.schemas import AnalysisResult, ReportArtifacts
from prostr_alfo.reporting.plots import plot_plddt_profile
from prostr_alfo.reporting.pymol import write_pymol_script
from prostr_alfo.reporting.report import generate_reports, load_references
from prostr_alfo.structure.alphafold import AlphaFoldClient, StructureNotFoundError
from prostr_alfo.structure.mutant import create_mutant_proxy_structure
from prostr_alfo.utils.paths import slugify
from prostr_alfo.variants.mutations import create_variant_registry, parse_mutations


@dataclass(slots=True)
class AnalysisRequest:
    uniprot_id: str | None = None
    sequence: str | None = None
    fasta_path: Path | None = None
    mutations: str | None = None
    label: str | None = None
    refresh_structure: bool = False


def run_analysis(request: AnalysisRequest, settings: Settings | None = None) -> AnalysisResult:
    """Run the full analysis workflow."""

    active_settings = settings or get_settings()
    active_settings.ensure_directories()

    protein_input = resolve_protein_input(
        uniprot_id=request.uniprot_id,
        sequence=request.sequence,
        fasta_path=request.fasta_path,
    )
    parsed_mutations = parse_mutations(request.mutations)

    label = request.label or protein_input.uniprot_id or protein_input.source_type
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_slug = slugify(f"{label}_{timestamp}")
    run_directory = active_settings.reports_dir / run_slug
    run_directory.mkdir(parents=True, exist_ok=True)

    variants_directory = active_settings.variants_dir / run_slug
    variants, registry_path = create_variant_registry(
        sequence=protein_input.sequence,
        mutations=parsed_mutations,
        output_dir=variants_directory,
        prefix=protein_input.uniprot_id or protein_input.source_type,
    )

    structure_analysis = None
    structure_available = False
    structure_note = None
    plot_path = None
    pymol_script_path = None
    mutant_structure_path = None

    if protein_input.uniprot_id:
        client = AlphaFoldClient(settings=active_settings)
        try:
            structure_path = client.get_structure_path(protein_input.uniprot_id, refresh=request.refresh_structure)
            structure_analysis = analyze_structure(structure_path, protein_input.sequence)
            structure_available = True
            mutant_structure_path = create_mutant_proxy_structure(
                structure_analysis=structure_analysis,
                mutations=parsed_mutations,
                output_path=run_directory / "mutant_proxy.pdb",
            )
            if mutant_structure_path is not None:
                structure_analysis.notes.append(
                    "A mutant proxy PDB was generated from WT coordinates without side-chain repacking or structural relaxation."
                )
            plot_path = plot_plddt_profile(
                structure_analysis=structure_analysis,
                mutations=parsed_mutations,
                output_path=run_directory / "plddt_profile.png",
            )
            pymol_script_path = write_pymol_script(
                structure_analysis=structure_analysis,
                mutations=parsed_mutations,
                mutant_structure_path=mutant_structure_path,
                output_path=run_directory / "view_structure.pml",
            )
        except StructureNotFoundError as exc:
            structure_note = str(exc)
        except Exception as exc:
            structure_note = f"Structure analysis failed: {exc}"
    else:
        structure_note = "No UniProt accession was supplied, so AlphaFold DB download was skipped."

    if structure_analysis is not None and structure_note:
        structure_analysis.notes.append(structure_note)

    mutation_assessments = assess_mutations(parsed_mutations, structure_analysis)
    variant_evidence = (
        collect_variant_evidence(
            accession=protein_input.uniprot_id,
            mutations=parsed_mutations,
            settings=active_settings,
        )
        if protein_input.uniprot_id and parsed_mutations
        else []
    )
    references = load_references(active_settings.references_path)
    markdown_path, html_path = generate_reports(
        protein_input=protein_input,
        structure_analysis=structure_analysis,
        variants=variants,
        mutation_assessments=mutation_assessments,
        variant_evidence=variant_evidence,
        references=references,
        plot_path=plot_path,
        pymol_script_path=pymol_script_path,
        registry_path=registry_path,
        output_dir=run_directory,
        structure_note=structure_note,
        mutant_structure_path=mutant_structure_path,
    )

    return AnalysisResult(
        protein_input=protein_input,
        structure_available=structure_available,
        structure_analysis=structure_analysis,
        variants=variants,
        mutation_assessments=mutation_assessments,
        variant_evidence=variant_evidence,
        report_artifacts=ReportArtifacts(
            markdown_path=markdown_path,
            html_path=html_path,
            plot_path=plot_path,
            pymol_script_path=pymol_script_path,
            registry_path=registry_path,
            mutant_structure_path=mutant_structure_path,
        ),
        references=references,
        run_directory=run_directory,
    )
