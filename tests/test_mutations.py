from pathlib import Path

import pytest

from prostr_alfo.variants.mutations import apply_mutations, create_variant_registry, parse_mutations


def test_parse_mutations_accepts_multiple_tokens() -> None:
    mutations = parse_mutations("A1V;C3S")
    assert [mutation.code for mutation in mutations] == ["A1V", "C3S"]


def test_parse_mutations_rejects_duplicate_positions() -> None:
    with pytest.raises(ValueError, match="provided more than once"):
        parse_mutations("A1V;A1G")


def test_apply_mutations_validates_wild_type_residue() -> None:
    mutations = parse_mutations("A1V")
    assert apply_mutations("ACDE", mutations) == "VCDE"

    mismatch = parse_mutations("C1V")
    with pytest.raises(ValueError, match="expects wild-type residue"):
        apply_mutations("ACDE", mismatch)


def test_create_variant_registry_writes_wt_and_cumulative_variants(tmp_path: Path) -> None:
    variants, registry_path = create_variant_registry(
        sequence="ACDEFG",
        mutations=parse_mutations("A1V;C2S"),
        output_dir=tmp_path,
        prefix="TEST",
    )

    assert [variant.name for variant in variants] == ["WT", "A1V", "A1V_C2S"]
    assert registry_path.exists()
    assert (tmp_path / "TEST_WT.fasta").exists()
    assert (tmp_path / "TEST_A1V.fasta").exists()
    assert (tmp_path / "TEST_A1V_C2S.fasta").exists()
