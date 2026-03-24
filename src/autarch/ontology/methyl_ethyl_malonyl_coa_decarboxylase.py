"""methyl/ethyl malonyl-CoA decarboxylase activity.

Catalysis of the reaction: (S)-methylmalonyl-CoA + H+ = CO2 + propanoyl-CoA or (2S)-ethylmalonyl-CoA + H+ = butanoyl-CoA + CO2.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_METHYLMALONYL_COA = "CHEBI:57326"
CHEBI_ETHYLMALONYL_COA = "CHEBI:60909"
CHEBI_PROPIONYL_COA = "CHEBI:57392"
CHEBI_BUTYRYL_COA = "CHEBI:57371"
CHEBI_CO2 = "CHEBI:16526"

PAIRINGS = {
    CHEBI_METHYLMALONYL_COA: CHEBI_PROPIONYL_COA,
    CHEBI_ETHYLMALONYL_COA: CHEBI_BUTYRYL_COA,
}


class MethylEthylMalonylCoADecarboxylase(ReactionClass):
    """methyl/ethyl malonyl-CoA decarboxylase activity.

    Catalysis of the reaction: (S)-methylmalonyl-CoA + H+ = CO2 + propanoyl-CoA or (2S)-ethylmalonyl-CoA + H+ = butanoyl-CoA + CO2.
    """

    GO_ID = "GO:0004492"
    EC_NUMBER_PREFIX = "4.1.1.94"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        substrate = next((chebi for chebi in PAIRINGS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(is_member=False, explanation="No methyl- or ethylmalonyl-CoA substrate detected")

        expected_product = PAIRINGS[substrate]
        if expected_product not in right_chebis or CHEBI_CO2 not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing the expected acyl-CoA product or carbon dioxide from methyl/ethylmalonyl-CoA decarboxylation",
            )

        substantive_left = {chebi for chebi in left_chebis if chebi not in {CHEBI_H_PLUS}}
        substantive_right = {chebi for chebi in right_chebis if chebi not in {CHEBI_H_PLUS}}
        if substantive_left != {substrate} or substantive_right != {expected_product, CHEBI_CO2}:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside methyl/ethylmalonyl-CoA decarboxylase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Methyl/ethyl malonyl-CoA decarboxylase: decarboxylation of a branched malonyl-CoA to the corresponding short-chain acyl-CoA",
        )
