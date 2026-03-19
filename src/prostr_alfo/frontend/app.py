"""Streamlit web interface for prostr-alfo."""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from prostr_alfo.analysis.pipeline import AnalysisRequest, run_analysis
from prostr_alfo.io.fasta import parse_fasta_text
from prostr_alfo.io.input import infer_uniprot_id_from_header

try:
    import py3Dmol
except ImportError:  # pragma: no cover - frontend fallback
    py3Dmol = None


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


def render_structure_panel(
    *,
    structure_path: Path | None,
    mutation_positions: list[str],
    low_confidence: list[str],
    disulfide_positions: list[str],
    mutant_mode: bool,
) -> None:
    """Render one interactive structure panel."""

    if py3Dmol is None:
        st.warning("py3Dmol is not installed, so the interactive 3D panel is disabled.")
        return

    if structure_path is None or not structure_path.exists():
        st.info("This structure view is not available for the current run.")
        return

    pdb_text = structure_path.read_text(encoding="utf-8")
    view = py3Dmol.view(width=920, height=560)
    view.addModel(pdb_text, "pdb")
    if mutant_mode:
        view.setStyle({}, {"cartoon": {"color": "lightblue"}})
    else:
        view.setStyle({}, {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 0, "max": 100}}})

    if low_confidence:
        view.addStyle({"resi": ",".join(low_confidence)}, {"stick": {"color": "orange", "radius": 0.18}})

    if mutation_positions:
        mutation_color = "red" if mutant_mode else "magenta"
        view.addStyle({"resi": ",".join(mutation_positions)}, {"sphere": {"color": mutation_color, "radius": 0.75}})

    if disulfide_positions:
        view.addStyle({"resi": ",".join(disulfide_positions)}, {"stick": {"color": "yellow", "radius": 0.25}})

    view.zoomTo()
    components.html(view._make_html(), height=580)


def collect_request() -> tuple[AnalysisRequest, bool]:
    """Collect one analysis request from the UI."""

    with st.form("analysis_form"):
        input_mode = st.radio("Input type", ["UniProt ID", "Raw sequence", "FASTA upload"], horizontal=True)
        uniprot_id = st.text_input("UniProt ID", placeholder="P05067").strip().upper()
        mutation_text = st.text_input("Mutations", placeholder="A23V;C105S").strip()
        label = st.text_input("Run label", placeholder="optional-analysis-name").strip() or None

        if input_mode == "UniProt ID":
            request = AnalysisRequest(uniprot_id=uniprot_id or None, mutations=mutation_text or None, label=label)
            submitted = st.form_submit_button("Run Analysis", type="primary", use_container_width=True)
            return request, submitted

        if input_mode == "Raw sequence":
            sequence = st.text_area("Protein sequence", height=180, placeholder="MSTNPKPQR...")
            request = AnalysisRequest(
                uniprot_id=uniprot_id or None,
                sequence=sequence or None,
                mutations=mutation_text or None,
                label=label,
            )
            submitted = st.form_submit_button("Run Analysis", type="primary", use_container_width=True)
            return request, submitted

        uploaded = st.file_uploader("Upload FASTA", type=["fasta", "fa", "faa", "txt"])
        sequence = None
        inferred_uniprot = None
        if uploaded is not None:
            header, sequence = parse_fasta_text(uploaded.getvalue().decode("utf-8"))
            inferred_uniprot = infer_uniprot_id_from_header(header)
            st.caption(f"Loaded FASTA header: {header}")
            if inferred_uniprot and not uniprot_id:
                st.caption(f"Inferred UniProt accession from FASTA header: {inferred_uniprot}")
        request = AnalysisRequest(
            uniprot_id=uniprot_id or inferred_uniprot,
            sequence=sequence,
            mutations=mutation_text or None,
            label=label,
        )
        submitted = st.form_submit_button("Run Analysis", type="primary", use_container_width=True)
        return request, submitted


def render_result(result) -> None:
    """Render analysis outputs."""

    variant_evidence = getattr(result, "variant_evidence", [])
    mutant_structure_path = getattr(result.report_artifacts, "mutant_structure_path", None)
    submitted_mutations = ", ".join(assessment.mutation.code for assessment in result.mutation_assessments) or "none"

    st.success(f"Analysis finished in {result.run_directory}")
    st.caption(f"Submitted mutations captured by the backend: {submitted_mutations}")

    if result.structure_analysis is not None:
        metric_a, metric_b, metric_c, metric_d = st.columns(4)
        metric_a.metric("Mean pLDDT", result.structure_analysis.mean_plddt)
        metric_b.metric("pLDDT > 70", f"{result.structure_analysis.percent_above_70}%")
        metric_c.metric("pLDDT > 90", f"{result.structure_analysis.percent_above_90}%")
        metric_d.metric("Surface proxy", f"{result.structure_analysis.surface_fraction}%")
    else:
        st.warning("Structure-based statistics are unavailable for this run.")

    tab_overview, tab_mutations, tab_evidence, tab_structure, tab_report = st.tabs(
        ["Overview", "Mutations", "Variant Evidence", "3D Structures", "Report"]
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

    with tab_evidence:
        if variant_evidence:
            for evidence in variant_evidence:
                st.markdown(f"### {evidence.mutation.code}")
                info_a, info_b, info_c = st.columns(3)
                info_a.metric("Exact curated match", "Yes" if evidence.matched else "No")
                info_b.metric("Uncertainty", evidence.uncertainty)
                info_c.metric("Clinical significance", ", ".join(evidence.clinical_significances) or "Unavailable")
                st.write(evidence.evidence_summary)
                st.write(f"Phenotype sufficiency note: {evidence.phenotype_sufficiency or 'Unavailable'}")
                if evidence.frequency_summary:
                    st.write(f"Frequency summary: {evidence.frequency_summary}")
                if evidence.pathways:
                    st.markdown("**Potentially affected pathways**")
                    st.dataframe(
                        [
                            {"Pathway": pathway.name, "Source": pathway.source, "ID": pathway.pathway_id, "URL": pathway.url}
                            for pathway in evidence.pathways
                        ],
                        use_container_width=True,
                    )
                if evidence.diseases:
                    st.markdown("**Disease associations**")
                    st.dataframe(
                        [
                            {
                                "Disease": disease.name,
                                "Description": disease.description or "No description available",
                                "PMIDs": ", ".join(disease.evidence_pmids[:5]),
                            }
                            for disease in evidence.diseases
                        ],
                        use_container_width=True,
                    )
                if evidence.literature:
                    st.markdown("**Linked papers**")
                    for paper in evidence.literature:
                        st.markdown(
                            f"- [{paper.title}]({paper.url}) ({paper.journal or 'Journal unavailable'}, {paper.publication_date or 'date unavailable'})"
                        )
                st.divider()
        else:
            st.info("No external variant evidence was available for this run.")

    with tab_structure:
        mutation_positions = [str(assessment.mutation.position) for assessment in result.mutation_assessments]
        low_confidence = [
            str(residue.sequence_position)
            for residue in (result.structure_analysis.residues if result.structure_analysis is not None else [])
            if residue.plddt is not None and residue.plddt < 50
        ]
        disulfide_positions = [
            str(position)
            for bond in (result.structure_analysis.disulfide_bonds if result.structure_analysis is not None else [])
            for position in (bond.residue_a, bond.residue_b)
        ]
        column_wt, column_mutant = st.columns(2)
        with column_wt:
            st.markdown("#### Wild-type structure")
            render_structure_panel(
                structure_path=result.structure_analysis.structure_path if result.structure_analysis is not None else None,
                mutation_positions=mutation_positions,
                low_confidence=low_confidence,
                disulfide_positions=disulfide_positions,
                mutant_mode=False,
            )
        with column_mutant:
            st.markdown("#### Mutant proxy structure")
            st.caption("WT coordinates with residue identity remapped at the mutation site. No side-chain repacking is performed.")
            render_structure_panel(
                structure_path=mutant_structure_path,
                mutation_positions=mutation_positions,
                low_confidence=[],
                disulfide_positions=[],
                mutant_mode=True,
            )

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
    if mutant_structure_path is not None:
        st.download_button(
            "Download mutant proxy PDB",
            data=mutant_structure_path.read_bytes(),
            file_name=mutant_structure_path.name,
            mime="chemical/x-pdb",
        )

    with st.expander("JSON summary"):
        st.code(json.dumps(result.to_dict(), indent=2), language="json")


st.subheader("Run analysis")
request, submitted = collect_request()
if submitted:
    try:
        if not any([request.uniprot_id, request.sequence, request.fasta_path]):
            raise ValueError("Please provide a UniProt ID, a raw sequence, or a FASTA upload.")
        with st.spinner("Running analysis..."):
            analysis_result = run_analysis(request)
        render_result(analysis_result)
    except Exception as exc:
        st.error(str(exc))
