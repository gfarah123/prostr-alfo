from prostr_alfo.structure.mapping import build_sequence_to_structure_map


def test_mapping_handles_missing_internal_residue() -> None:
    mapping = build_sequence_to_structure_map("ACDEFG", "ACEFG")
    assert mapping == {1: 1, 2: 2, 3: None, 4: 3, 5: 4, 6: 5}


def test_mapping_handles_missing_n_terminal_residue() -> None:
    mapping = build_sequence_to_structure_map("MACDE", "ACDE")
    assert mapping == {1: None, 2: 1, 3: 2, 4: 3, 5: 4}
