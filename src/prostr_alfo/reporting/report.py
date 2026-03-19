"""Markdown and HTML report generation."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import markdown

from prostr_alfo.models.schemas import MutationAssessment, ProteinInput, StructureAnalysis, VariantRecord


def load_references(path: Path) -> list[dict]:
    """Load structured references from disk."""

    return json.loads(path.read_text(encoding="utf-8"))


def _format_structure_summary(structure_analysis: StructureAnalysis | None, structure_note: str | None) -> str:
    if structure_analysis is None:
        note = structure_note or "No structure analysis was available."
        return f"- Structure status: unavailable\n- Note: {note}\n"

    low_regions = ", ".join(
        f"{region.start}-{region.end} (mean {region.mean_plddt:.1f})"
        for region in structure_analysis.low_confidence_regions
    ) or "None detected"
    disulfides = ", ".join(
        f"Cys {bond.residue_a}-Cys {bond.residue_b} ({bond.distance:.2f} A)"
        for bond in structure_analysis.disulfide_bonds
    ) or "No putative disulfides detected"
    secondary = ", ".join(f"{code}: {count}" for code, count in sorted(structure_analysis.secondary_structure_counts.items())) or "Unavailable"

    lines = [
        f"- Chains: {', '.join(structure_analysis.chain_ids)}",
        f"- Mean pLDDT: {structure_analysis.mean_plddt}",
        f"- Residues with pLDDT > 70: {structure_analysis.percent_above_70}%",
        f"- Residues with pLDDT > 90: {structure_analysis.percent_above_90}%",
        f"- Surface proxy fraction: {structure_analysis.surface_fraction}%",
        f"- Low-confidence regions: {low_regions}",
        f"- Secondary structure summary: {secondary}",
        f"- Disulfide heuristic: {disulfides}",
    ]
    if structure_analysis.notes:
        lines.append(f"- Notes: {'; '.join(structure_analysis.notes)}")
    return "\n".join(lines) + "\n"


def _format_mutation_table(assessments: list[MutationAssessment]) -> str:
    if not assessments:
        return "No mutations were requested.\n"

    lines = [
        "| Mutation | Mapped | pLDDT | Burial | Contacts | Charge delta | Hydrophobicity delta | Score | Level |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for assessment in assessments:
        lines.append(
            "| {mutation} | {mapped} | {plddt} | {burial} | {contacts} | {charge} | {hydro} | {score} | {level} |".format(
                mutation=assessment.mutation.code,
                mapped="Yes" if assessment.mapped else "No",
                plddt=f"{assessment.residue_plddt:.1f}" if assessment.residue_plddt is not None else "-",
                burial=f"{assessment.burial_score:.2f}" if assessment.burial_score is not None else "-",
                contacts=assessment.local_contact_count if assessment.local_contact_count is not None else "-",
                charge=f"{assessment.charge_delta:+.1f}" if assessment.charge_delta is not None else "-",
                hydro=f"{assessment.hydrophobicity_delta:+.2f}" if assessment.hydrophobicity_delta is not None else "-",
                score=f"{assessment.impact_score:.2f}" if assessment.impact_score is not None else "-",
                level=assessment.impact_level or "-",
            )
        )
    return "\n".join(lines) + "\n"


def _format_variant_table(variants: list[VariantRecord]) -> str:
    lines = [
        "| Variant | Mutations | FASTA |",
        "| --- | --- | --- |",
    ]
    for variant in variants:
        mutation_text = ", ".join(mutation.code for mutation in variant.mutations) or "WT"
        lines.append(f"| {variant.name} | {mutation_text} | `{variant.fasta_path.name}` |")
    return "\n".join(lines) + "\n"


def _format_references(references: list[dict]) -> str:
    lines = []
    for reference in references:
        citation = (
            f"- **{reference['id']}**: {reference['authors']} ({reference['year']}). "
            f"{reference['title']}. *{reference['journal']}*. "
            f"[DOI: {reference['doi']}]({reference['url']})"
        )
        lines.append(citation)
    return "\n".join(lines) + "\n"


def _build_markdown_report(
    *,
    protein_input: ProteinInput,
    structure_analysis: StructureAnalysis | None,
    variants: list[VariantRecord],
    mutation_assessments: list[MutationAssessment],
    references: list[dict],
    plot_path: Path | None,
    pymol_script_path: Path | None,
    registry_path: Path,
    structure_note: str | None,
) -> str:
    plot_line = f"![pLDDT profile]({plot_path.name})\n" if plot_path is not None else "Plot unavailable.\n"
    pymol_line = pymol_script_path.name if pymol_script_path is not None else "Unavailable"

    mutation_interpretations = "\n".join(
        f"- **{assessment.mutation.code}**: {assessment.interpretation or assessment.message}"
        for assessment in mutation_assessments
    ) or "- No mutation interpretations available."

    return f"""# prostr-alfo Analysis Report

## Protein Summary

- Input source: {protein_input.source_type}
- Source value: {protein_input.source_value}
- UniProt accession: {protein_input.uniprot_id or "Not supplied"}
- Sequence length: {len(protein_input.sequence)}

## Structure and Quality Features

{_format_structure_summary(structure_analysis, structure_note)}

## pLDDT Plot

{plot_line}

## Variant Registry

{_format_variant_table(variants)}

Registry JSON: `{registry_path.name}`

## Mutation Analysis

{_format_mutation_table(mutation_assessments)}

### Interpretation

{mutation_interpretations}

## Methods

- Sequence input was normalized from UniProt FASTA, uploaded FASTA, or raw amino-acid text.
- When a UniProt accession was available, AlphaFold DB metadata was queried and a PDB file was cached locally.
- Structure parsing used Biopython. pLDDT values were extracted from atom B-factors in the AlphaFold PDB.
- Local packing was approximated using C-alpha neighbor counts within 10 A, and surface exposure was approximated from local density.
- Mutation impact scoring combines residue confidence, burial proxy, local contact density, side-chain volume change, charge change, hydrophobicity change, and special-residue rules.
- Variant FASTA files were produced for WT and cumulative mutation combinations.
- PyMOL visualization script: `{pymol_line}`

## Limitations

- The mutation impact score is a heuristic ranking aid and is not a thermodynamic free-energy predictor.
- Sequence-only analyses without a resolvable UniProt accession cannot retrieve AlphaFold structures automatically.
- AlphaFold structures are predictions, and flexible or disordered regions can still be unreliable even when coordinates are present.
- DSSP annotations require an external DSSP installation; when unavailable, the report omits secondary structure labels.
- Contact and surface metrics use coarse geometric proxies rather than full solvent-accessible surface area or energetic calculations.

## References

{_format_references(references)}
"""


def _build_html(markdown_text: str, title: str) -> str:
    body = markdown.markdown(markdown_text, extensions=["tables", "fenced_code"])
    escaped_title = escape(title)
    return f"""<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{escaped_title}</title>
    <style>
      :root {{
        color-scheme: light;
        --bg: #f5f7fb;
        --panel: #ffffff;
        --ink: #0f172a;
        --muted: #475569;
        --line: #dbe3f0;
        --accent: #0f766e;
      }}
      body {{
        background: radial-gradient(circle at top, #e0f2fe 0%, var(--bg) 45%, #eef2ff 100%);
        color: var(--ink);
        font-family: "Segoe UI", "Helvetica Neue", sans-serif;
        margin: 0;
        padding: 2rem 1rem;
      }}
      main {{
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 18px;
        box-shadow: 0 18px 45px rgba(15, 23, 42, 0.08);
        margin: 0 auto;
        max-width: 1080px;
        padding: 2rem;
      }}
      h1, h2, h3 {{
        color: #0f172a;
      }}
      p, li {{
        color: var(--muted);
        line-height: 1.6;
      }}
      a {{
        color: var(--accent);
      }}
      table {{
        border-collapse: collapse;
        margin: 1rem 0;
        width: 100%;
      }}
      th, td {{
        border: 1px solid var(--line);
        padding: 0.65rem;
        text-align: left;
      }}
      th {{
        background: #eff6ff;
      }}
      img {{
        border-radius: 12px;
        max-width: 100%;
      }}
      code {{
        background: #eff6ff;
        border-radius: 6px;
        padding: 0.12rem 0.32rem;
      }}
    </style>
  </head>
  <body>
    <main>{body}</main>
  </body>
</html>
"""


def generate_reports(
    *,
    protein_input: ProteinInput,
    structure_analysis: StructureAnalysis | None,
    variants: list[VariantRecord],
    mutation_assessments: list[MutationAssessment],
    references: list[dict],
    plot_path: Path | None,
    pymol_script_path: Path | None,
    registry_path: Path,
    output_dir: Path,
    structure_note: str | None,
) -> tuple[Path, Path]:
    """Generate markdown and HTML reports in one output directory."""

    markdown_text = _build_markdown_report(
        protein_input=protein_input,
        structure_analysis=structure_analysis,
        variants=variants,
        mutation_assessments=mutation_assessments,
        references=references,
        plot_path=plot_path,
        pymol_script_path=pymol_script_path,
        registry_path=registry_path,
        structure_note=structure_note,
    )
    markdown_path = output_dir / "report.md"
    html_path = output_dir / "report.html"
    markdown_path.write_text(markdown_text, encoding="utf-8")
    html_path.write_text(_build_html(markdown_text, "prostr-alfo report"), encoding="utf-8")
    return markdown_path, html_path
