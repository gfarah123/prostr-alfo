from pathlib import Path

from prostr_alfo.models.schemas import StructureAnalysis
from prostr_alfo.structure.mutant import create_mutant_proxy_structure
from prostr_alfo.variants.mutations import parse_mutations


def test_create_mutant_proxy_structure_renames_residue_and_prunes_atoms(tmp_path: Path) -> None:
    wt_pdb = tmp_path / "wt.pdb"
    wt_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1      11.104  13.207   9.147  1.00 90.00           N",
                "ATOM      2  CA  ALA A   1      12.560  13.207   9.147  1.00 90.00           C",
                "ATOM      3  C   ALA A   1      13.000  14.620   9.147  1.00 90.00           C",
                "ATOM      4  O   ALA A   1      12.400  15.660   9.147  1.00 90.00           O",
                "ATOM      5  CB  ALA A   1      13.100  12.100   8.200  1.00 90.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    structure_analysis = StructureAnalysis(
        structure_path=wt_pdb,
        chain_ids=["A"],
        structure_sequence="A",
        residues=[],
        mean_plddt=90.0,
        percent_above_70=100.0,
        percent_above_90=100.0,
        sequence_to_structure={1: 1},
    )

    output_path = tmp_path / "mutant.pdb"
    mutated = create_mutant_proxy_structure(
        structure_analysis=structure_analysis,
        mutations=parse_mutations("A1G"),
        output_path=output_path,
    )

    assert mutated == output_path
    mutated_text = output_path.read_text(encoding="utf-8")
    assert "GLY" in mutated_text
    assert " CB " not in mutated_text
