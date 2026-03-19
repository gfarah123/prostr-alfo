# prostr-alfo

`prostr-alfo` is a Python and Streamlit application for protein structure and mutation analysis built around AlphaFold DB. It accepts UniProt accessions, raw protein sequences, or FASTA input, generates mutation variants, extracts confidence and geometric structure features, scores point mutations with transparent heuristics, and produces reproducible Markdown and HTML reports.

## Features

- UniProt, raw sequence, and FASTA input support
- AlphaFold DB download and local PDB caching for UniProt-backed analyses
- pLDDT extraction from AlphaFold structures
- Structure feature extraction:
  - mean pLDDT
  - percent of residues above 70 and 90
  - low-confidence regions below 50
  - optional DSSP secondary structure assignment
  - Cys-Cys disulfide distance heuristic
  - CA-density surface and burial proxy
- Mutation engine with validation and FASTA output
- Cumulative variant registry generation:
  - `WT`
  - `A23V`
  - `A23V_C105S`
- Mutation assessment using WT structure context
- Heuristic mutation impact score with interpretable components
- Curated external evidence lookup for pathways, disease associations, clinical significance, and linked papers
- Report generation in Markdown and HTML
- PyMOL script generation
- Interactive Streamlit interface with WT and mutant proxy 3D views
- Pytest coverage for parsing, mutation application, caching, and mapping edge cases

## Installation

1. Create and activate a virtual environment.
2. Install the project in editable mode.

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install pytest
```

You can also install from `requirements.txt`:

```powershell
python -m pip install -r requirements.txt
```

## CLI Usage

Run the main workflow with Typer:

```powershell
prostr-alfo analyze --uniprot P69905 --mutations "M1V"
```

```powershell
prostr-alfo analyze --sequence "ACDEFGHIK" --mutations "A1V"
```

```powershell
prostr-alfo analyze --fasta .\example.fasta --mutations "A23V;C105S"
```

Useful options:

- `--label`: custom output label
- `--refresh-structure`: force a fresh AlphaFold download
- `--json-output`: print the JSON payload to the terminal

## Frontend Usage

Run the Streamlit app:

```powershell
streamlit run src/prostr_alfo/frontend/app.py
```

The frontend supports:

- UniProt ID entry
- raw sequence input
- FASTA upload
- mutation input
- pLDDT plot display
- mutation summary table
- variant evidence tab with pathways, disease associations, frequencies, and linked papers
- side-by-side WT and mutant proxy structure views with mutation highlighting
- report preview
- downloads for HTML report, FASTA, PyMOL script, and variant registry
- mutant proxy PDB download when available

## Example Workflow

1. Start with a UniProt accession such as `P69905`.
2. Run:

```powershell
prostr-alfo analyze --uniprot P69905 --mutations "M1V"
```

3. Review outputs in `data/reports/<run_id>/`:
   - `report.md`
   - `report.html`
   - `plddt_profile.png`
   - `mutant_proxy.pdb`
   - `view_structure.pml`
4. Review generated FASTA files in `data/variants/<run_id>/`.
5. Open the HTML report in a browser or use the Streamlit app for interactive inspection.

## Project Structure

```text
prostr-alfo/
  README.md
  pyproject.toml
  requirements.txt
  .gitignore
  src/
    prostr_alfo/
      cli.py
      config.py
      models/
        schemas.py
      io/
        fasta.py
        input.py
      structure/
        alphafold.py
        mapping.py
        parser.py
      variants/
        mutations.py
      analysis/
        features.py
        mutation.py
        pipeline.py
        scoring.py
      reporting/
        plots.py
        pymol.py
        references.json
        report.py
      frontend/
        app.py
      utils/
        biology.py
        paths.py
  data/
    cache/
    reports/
    structures/
    variants/
  tests/
```

## Output Layout

- `data/structures/`: cached AlphaFold PDB files
- `data/variants/`: generated FASTA files and variant registry JSON
- `data/reports/`: Markdown report, HTML report, plot, and PyMOL script

## Methods Summary

- UniProt sequences are retrieved from the UniProt REST FASTA endpoint.
- AlphaFold structures are retrieved from AlphaFold DB metadata and stored locally.
- pLDDT is read from atom B-factors in AlphaFold PDB files.
- Low-confidence regions are contiguous segments with pLDDT below 50.
- Surface exposure is approximated by CA neighbor density within 10 A.
- Mutation impact scoring combines:
  - local pLDDT
  - burial proxy
  - local contact density
  - side-chain volume change proxy
  - charge change
  - hydrophobicity change
  - special residue rules for Gly, Pro, and Cys
- Curated variant evidence uses:
  - Reactome for protein-level pathway context
  - EBI Proteins variation records for exact amino-acid substitutions when available
  - PubMed citation metadata for linked papers

## Scientific References

The structured references are stored in [`src/prostr_alfo/reporting/references.json`](/Users/41768/Desktop/prostr-alfo/src/prostr_alfo/reporting/references.json).

- Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583-589 (2021). DOI: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- Tunyasuvunakool K, Adler J, Wu Z, et al. Highly accurate protein structure prediction for the human proteome. Nature 596, 590-596 (2021). DOI: [10.1038/s41586-021-03828-1](https://doi.org/10.1038/s41586-021-03828-1)
- Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers 22(12), 2577-2637 (1983). DOI: [10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)
- Stefl S, Nishi H, Petukh M, Panchenko AR, Alexov E. Molecular mechanisms of disease-causing missense mutations. Journal of Molecular Biology 425(21), 3919-3936 (2013). DOI: [10.1016/j.jmb.2013.02.041](https://doi.org/10.1016/j.jmb.2013.02.041)

## Limitations

- Full structure retrieval depends on having a resolvable UniProt accession. Raw sequence and FASTA inputs without a matching accession can still be analyzed, but structure-dependent steps are skipped.
- The mutation impact score is heuristic and should not be interpreted as a calibrated stability or pathogenicity predictor.
- The mutant structure is a WT-coordinate proxy without side-chain repacking, minimization, or a new structure prediction.
- The surface and contact metrics are geometric proxies, not molecular mechanics calculations.
- DSSP annotations require a local DSSP executable installation.
- Multi-chain biological assemblies are only handled at the single downloaded PDB level and are not reassembled from external complex predictions.
- Linked papers can be disease- or gene-level support rather than dedicated residue-specific functional studies, and the report marks uncertainty explicitly when evidence is indirect or conflicting.

## Future Improvements

- Add explicit solvent-accessible surface area calculations
- Add side-chain repacking or external stability engines
- Support AlphaFold multimer and complex workflows
- Add richer residue environment plots and contact maps
- Add exportable CSV summaries
- Add background job handling for larger proteins in the web app

## Git Workflow

This project should not be developed directly on `main`.

### First-time repository setup

```powershell
git init
git remote add origin https://github.com/gfarah123/prostr-alfo.git
git checkout -b feature/initial-setup
git add .
git commit -m "Initial project setup"
git push -u origin feature/initial-setup
git checkout main
git merge feature/initial-setup
git push origin main
```

### Ongoing work

```powershell
git checkout main
git pull origin main
git checkout -b feature/<name>
git add .
git commit -m "Describe the change"
git checkout main
git merge feature/<name>
git push origin main
```

## Testing

Run the test suite with:

```powershell
python -m pytest
```

## Windows Notes

- All project paths are handled with `pathlib` for Windows compatibility.
- The Streamlit app can be launched directly with `streamlit run src/prostr_alfo/frontend/app.py`.
- Cached structures and generated reports are stored inside the local project folder.
