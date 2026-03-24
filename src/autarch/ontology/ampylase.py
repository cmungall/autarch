"""AMPylase activity.

Catalysis of the reaction: ATP + protein = diphosphate + adenylyl-protein; mediates the addition of an adenylyl (adenosine 5'-monophosphate; AMP group) to L-serine, L-threonine, and L-tyrosine residues in target proteins.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ATP, CHEBI_DIPHOSPHATE
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

AMPYLATION_PAIRS = {
    "CHEBI:29999": "CHEBI:142516",
    "CHEBI:30013": "CHEBI:138113",
    "CHEBI:46858": "CHEBI:83624",
}


class AMPylase(ReactionClass):
    """AMPylase activity.

    Catalysis of the reaction: ATP + protein = diphosphate + adenylyl-protein; mediates the addition of an adenylyl (adenosine 5'-monophosphate; AMP group) to L-serine, L-threonine, and L-tyrosine residues in target proteins.
    """

    GO_ID = "GO:0070733"
    EC_NUMBER_PREFIX = "2.7.7.108"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}

        if CHEBI_ATP not in left_chebis or CHEBI_DIPHOSPHATE not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Requires ATP substrate and diphosphate product",
            )

        substrate = next((chebi for chebi in AMPYLATION_PAIRS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(
                is_member=False,
                explanation="No protein-bound serine, threonine, or tyrosine AMPylation substrate detected",
            )
        product = AMPYLATION_PAIRS[substrate]
        if product not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing the corresponding adenylylated protein product",
            )

        substantive_left = left_chebis - {CHEBI_ATP, substrate}
        substantive_right = right_chebis - {CHEBI_DIPHOSPHATE, product}
        if substantive_left or substantive_right:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside AMPylase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="AMPylase: ATP-dependent transfer of AMP onto a protein-bound serine, threonine, or tyrosine residue",
        )
