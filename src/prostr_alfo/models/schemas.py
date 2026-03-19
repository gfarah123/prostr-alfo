"""Typed dataclasses shared across modules."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class ProteinInput:
    source_type: str
    source_value: str
    sequence: str
    uniprot_id: str | None = None
    fasta_header: str | None = None


@dataclass(slots=True, frozen=True)
class Mutation:
    wild_type: str
    position: int
    mutant: str

    @property
    def code(self) -> str:
        return f"{self.wild_type}{self.position}{self.mutant}"


@dataclass(slots=True)
class ResidueFeature:
    sequence_position: int
    structure_index: int | None
    chain_id: str
    residue_name: str
    plddt: float | None
    secondary_structure: str | None
    contact_density: int
    burial_score: float | None
    is_surface: bool
    x: float | None
    y: float | None
    z: float | None


@dataclass(slots=True)
class LowConfidenceRegion:
    start: int
    end: int
    mean_plddt: float


@dataclass(slots=True)
class DisulfideBond:
    residue_a: int
    residue_b: int
    distance: float


@dataclass(slots=True)
class StructureAnalysis:
    structure_path: Path
    chain_ids: list[str]
    structure_sequence: str
    residues: list[ResidueFeature]
    mean_plddt: float | None
    percent_above_70: float | None
    percent_above_90: float | None
    low_confidence_regions: list[LowConfidenceRegion] = field(default_factory=list)
    disulfide_bonds: list[DisulfideBond] = field(default_factory=list)
    secondary_structure_counts: dict[str, int] = field(default_factory=dict)
    surface_fraction: float | None = None
    notes: list[str] = field(default_factory=list)
    sequence_to_structure: dict[int, int | None] = field(default_factory=dict)


@dataclass(slots=True)
class MutationAssessment:
    mutation: Mutation
    mapped: bool
    message: str
    residue_plddt: float | None = None
    burial_score: float | None = None
    local_contact_count: int | None = None
    contact_change_proxy: float | None = None
    hydrophobicity_delta: float | None = None
    charge_delta: float | None = None
    special_residue_flags: list[str] = field(default_factory=list)
    neighboring_residues: list[str] = field(default_factory=list)
    impact_score: float | None = None
    impact_level: str | None = None
    interpretation: str | None = None


@dataclass(slots=True)
class VariantRecord:
    name: str
    mutations: list[Mutation]
    sequence: str
    fasta_path: Path


@dataclass(slots=True)
class ReportArtifacts:
    markdown_path: Path
    html_path: Path
    plot_path: Path | None
    pymol_script_path: Path | None
    registry_path: Path


@dataclass(slots=True)
class AnalysisResult:
    protein_input: ProteinInput
    structure_available: bool
    structure_analysis: StructureAnalysis | None
    variants: list[VariantRecord]
    mutation_assessments: list[MutationAssessment]
    report_artifacts: ReportArtifacts
    references: list[dict[str, Any]]
    run_directory: Path

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["run_directory"] = str(self.run_directory)
        payload["report_artifacts"] = {
            "markdown_path": str(self.report_artifacts.markdown_path),
            "html_path": str(self.report_artifacts.html_path),
            "plot_path": str(self.report_artifacts.plot_path) if self.report_artifacts.plot_path else None,
            "pymol_script_path": str(self.report_artifacts.pymol_script_path) if self.report_artifacts.pymol_script_path else None,
            "registry_path": str(self.report_artifacts.registry_path),
        }
        if self.structure_analysis is not None:
            payload["structure_analysis"]["structure_path"] = str(self.structure_analysis.structure_path)
        for variant in payload["variants"]:
            variant["fasta_path"] = str(variant["fasta_path"])
        return payload
