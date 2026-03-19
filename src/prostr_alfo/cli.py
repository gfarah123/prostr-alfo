"""Typer CLI entry points."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import typer

from prostr_alfo.analysis.pipeline import AnalysisRequest, run_analysis

app = typer.Typer(help="Protein structure and mutation analysis using AlphaFold DB.")


@app.command()
def analyze(
    uniprot: str | None = typer.Option(None, "--uniprot", help="UniProt accession."),
    sequence: str | None = typer.Option(None, "--sequence", help="Raw amino-acid sequence."),
    fasta: Path | None = typer.Option(None, "--fasta", exists=True, dir_okay=False, file_okay=True, help="FASTA file path."),
    mutations: str | None = typer.Option(None, "--mutations", help="Mutation expression such as A23V;C105S."),
    label: str | None = typer.Option(None, "--label", help="Optional run label."),
    refresh_structure: bool = typer.Option(False, "--refresh-structure", help="Force a fresh AlphaFold download."),
    json_output: bool = typer.Option(False, "--json-output", help="Print the full JSON result payload."),
) -> None:
    """Run the main analysis workflow."""

    request = AnalysisRequest(
        uniprot_id=uniprot,
        sequence=sequence,
        fasta_path=fasta,
        mutations=mutations,
        label=label,
        refresh_structure=refresh_structure,
    )
    result = run_analysis(request)

    typer.echo(f"Run directory: {result.run_directory}")
    typer.echo(f"Markdown report: {result.report_artifacts.markdown_path}")
    typer.echo(f"HTML report: {result.report_artifacts.html_path}")
    typer.echo(f"Structure available: {'yes' if result.structure_available else 'no'}")
    typer.echo(f"Variants generated: {len(result.variants)}")
    typer.echo(f"Mutation assessments: {len(result.mutation_assessments)}")

    if json_output:
        typer.echo(json.dumps(result.to_dict(), indent=2))


@app.command()
def web() -> None:
    """Launch the Streamlit frontend."""

    app_path = Path(__file__).resolve().parent / "frontend" / "app.py"
    command = [sys.executable, "-m", "streamlit", "run", str(app_path)]
    raise SystemExit(subprocess.call(command))


if __name__ == "__main__":
    app()
