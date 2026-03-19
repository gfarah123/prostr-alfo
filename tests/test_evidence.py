from prostr_alfo.analysis.evidence import match_variant_feature
from prostr_alfo.models.schemas import Mutation


def test_match_variant_feature_requires_exact_position_and_residue_change() -> None:
    dataset = {
        "features": [
            {
                "type": "VARIANT",
                "begin": "2",
                "wildType": "D",
                "mutatedType": "A",
            },
            {
                "type": "VARIANT",
                "begin": "2",
                "wildType": "D",
                "mutatedType": "G",
            },
        ]
    }

    matched = match_variant_feature(dataset, Mutation(wild_type="D", position=2, mutant="A"))
    assert matched is not None
    assert matched["mutatedType"] == "A"

    missing = match_variant_feature(dataset, Mutation(wild_type="D", position=3, mutant="A"))
    assert missing is None
