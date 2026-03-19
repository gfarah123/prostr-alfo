"""Streamlit web interface for prostr-alfo."""

from __future__ import annotations

import json
from pathlib import Path

import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from prostr_alfo.analysis.pipeline import AnalysisRequest, run_analysis
from prostr_alfo.io.fasta import parse_fasta_text
from prostr_alfo.io.input import infer_uniprot_id_from_header


st.set_page_config(page_title="prostr-alfo", layout="wide")

st.markdown(
    """
    <style>
      .stApp {
        background:
          radial-gradient(circle at top left, rgba(191, 219, 254, 0.85), transparent 26%),
          radial-gradient(circle at top right, rgba(253, 230, 138, 0.75), transparent 22%),
          linear-gradient(180deg, #f8fafc 0%, #eef2ff 100%);
      }
      .hero {
        background: rgba(255, 255, 255, 0.85);
        border: 1px solid rgba(148, 163, 184, 0.35);
        border-radius: 24px;
        box-shadow: 0 24px 60px rgba(15, 23, 42, 0.08);
        padding: 1.8rem 2rem;
        margin-bottom: 1.5rem;
      }
      .hero h1 {
        color: #0f172a;
        font-family: Georgia, "Times New Roman", serif;
        font-size: 2.4rem;
        margin-bottom: 0.2rem;
      }
      .hero p {
        color: #475569;
        font-size: 1.02rem;
        margin: 0;
      }
    </style>
    <div class="hero">
      <h1>prostr-alfo</h1>
      <p>AlphaFold-based protein and mutation analysis with reproducible reports.</p>
    </div>
    """,
    unsafe_allow_html=True,
)


def render_structure_view(result) -> None:
    """Render an interactive 3D structure view."""

    if not result.structure_available or result.structure_analysis is None:
        st.info("No AlphaFold structure is available for interactive visualization.")
        return

    pdb_text = Path(result.structure_analysis.structure_path).read_text(encoding="utf-8")
    view = py3Dmol.view(width=920, height=560)
    view.addModel(pdb_text, "pdb")
    view.setStyle({}, {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 0, "max": 100}}})

    low_confidence = [str(residue.sequence_position) for residue in result.structure_analysis.residues if residue.plddt is not None and residue.plddt < 50]
    if low_confidence:
        view.addStyle({"resi": ",".join(low_confidence)}, {"stick": {"color": "orange", "radius": 0.18}})

    mutation_positions = [str(assessment.mutation.position) for assessment in result.mutation_assessments]
    if mutation_positions:
        view.addStyle({"resi": ",".join(mutation_positions)}, {"sphere": {"color": "magenta", "radius": 0.75}})

    disulfide_positions = [
        str(position)
        for bond in result.structure_analysis.disulfide_bonds
        for position in (bond.residue_a, bond.residue_b)
    ]
    if disulfide_positions:
        view.addStyle({"resi": ",".join(disulfide_positions)}, {"stick": {"color": "yellow", "radius": 0.25}})

    view.zoomTo()
    components.html(view._make_html(), height=580)


def collect_request() -> AnalysisRequest:
    """Collect one analysis request from the UI."""

    input_mode = st.radio("Input type", ["UniProt ID", "Raw sequence", "FASTA upload"], horizontal=True)
    uniprot_id = st.text_input("UniProt ID", placeholder="P05067").strip().upper()
    mutation_text = st.text_input("Mutations", placeholder="A23V;C105S").strip()
    label = st.text_input("Run label", placeholder="optional-analysis-name").strip() or None

    if input_mode == "UniProt ID":
        return AnalysisRequest(uniprot_id=uniprot_id or None, mutations=mutation_text or None, label=label)

    if input_mode == "Raw sequence":
        sequence = st.text_area("Protein sequence", height=180, placeholder="MSTNPKPQR...")
        return AnalysisRequest(
            uniprot_id=uniprot_id or None,
            sequence=sequence or None,
            mutations=mutation_text or None,
            label=label,
        )

    uploaded = st.file_uploader("Upload FASTA", type=["fasta", "fa", "faa", "txt"])
    sequence = None
    inferred_uniprot = None
    if uploaded is not None:
        header, sequence = parse_fasta_text(uploaded.getvalue().decode("utf-8"))
        inferred_uniprot = infer_uniprot_id_from_header(header)
        st.caption(f"Loaded FASTA header: {header}")
        if inferred_uniprot and not uniprot_id:
            st.caption(f"Inferred UniProt accession from FASTA header: {inferred_uniprot}")
    return AnalysisRequest(
        uniprot_id=uniprot_id or inferred_uniprot,
        sequence=sequence,
        mutations=mutation_text or None,
        label=label,
    )


def render_result(result) -> None:
    """Render analysis outputs."""

    st.success(f"Analysis finished in {result.run_directory}")

    if result.structure_analysis is not None:
        metric_a, metric_b, metric_c, metric_d = st.columns(4)
        metric_a.metric("Mean pLDDT", result.structure_analysis.mean_plddt)
        metric_b.metric("pLDDT > 70", f"{result.structure_analysis.percent_above_70}%")
        metric_c.metric("pLDDT > 90", f"{result.structure_analysis.percent_above_90}%")
        metric_d.metric("Surface proxy", f"{result.structure_analysis.surface_fraction}%")
    else:
        st.warning("Structure-based statistics are unavailable for this run.")

    tab_overview, tab_mutations, tab_structure, tab_report = st.tabs(
        ["Overview", "Mutations", "3D Structure", "Report"]
    )

    with tab_overview:
        if result.report_artifacts.plot_path is not None:
            st.image(str(result.report_artifacts.plot_path), caption="pLDDT profile")
        st.subheader("Variants")
        st.dataframe(
            [
                {
                    "Variant": variant.name,
                    "Mutations": ", ".join(mutation.code for mutation in variant.mutations) or "WT",
                    "FASTA": variant.fasta_path.name,
                }
                for variant in result.variants
            ],
            use_container_width=True,
        )

    with tab_mutations:
        if result.mutation_assessments:
            st.dataframe(
                [
                    {
                        "Mutation": assessment.mutation.code,
                        "Mapped": assessment.mapped,
                        "pLDDT": assessment.residue_plddt,
                        "Burial": assessment.burial_score,
                        "Contacts": assessment.local_contact_count,
                        "Score": assessment.impact_score,
                        "Level": assessment.impact_level,
                        "Interpretation": assessment.interpretation or assessment.message,
                    }
                    for assessment in result.mutation_assessments
                ],
                use_container_width=True,
            )
        else:
            st.info("No mutations were requested.")

    with tab_structure:
        render_structure_view(result)

    with tab_report:
        markdown_text = Path(result.report_artifacts.markdown_path).read_text(encoding="utf-8")
        st.markdown(markdown_text)

    html_bytes = Path(result.report_artifacts.html_path).read_bytes()
    st.download_button("Download HTML report", data=html_bytes, file_name=result.report_artifacts.html_path.name, mime="text/html")

    registry_bytes = Path(result.report_artifacts.registry_path).read_bytes()
    st.download_button(
        "Download variant registry",
        data=registry_bytes,
        file_name=result.report_artifacts.registry_path.name,
        mime="application/json",
    )

    if result.variants:
        first_variant = result.variants[-1]
        st.download_button(
            "Download latest variant FASTA",
            data=first_variant.fasta_path.read_bytes(),
            file_name=first_variant.fasta_path.name,
            mime="text/plain",
        )

    if result.report_artifacts.pymol_script_path is not None:
        st.download_button(
            "Download PyMOL script",
            data=result.report_artifacts.pymol_script_path.read_bytes(),
            file_name=result.report_artifacts.pymol_script_path.name,
            mime="text/plain",
        )

    with st.expander("JSON summary"):
        st.code(json.dumps(result.to_dict(), indent=2), language="json")


st.subheader("Run analysis")
request = collect_request()
if st.button("Run Analysis", type="primary", use_container_width=True):
    try:
        if not any([request.uniprot_id, request.sequence, request.fasta_path]):
            raise ValueError("Please provide a UniProt ID, a raw sequence, or a FASTA upload.")
        with st.spinner("Running analysis..."):
            analysis_result = run_analysis(request)
        render_result(analysis_result)
    except Exception as exc:
        st.error(str(exc))
