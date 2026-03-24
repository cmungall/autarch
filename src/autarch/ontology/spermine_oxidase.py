"""spermine oxidase activity.

Catalysis of the reaction: H2O + O2 + spermine = 3-aminopropanal + H2O2 + spermidine. Weak activity with N(1)-acetylspermine. The Arabidopsis thaliana enzyme converts norspermine to norspermidine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H2O2, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_AMINOPROPANAL = "CHEBI:133427"

POLYAMINE_PRODUCT_PAIRS = {
    "CHEBI:45725": "CHEBI:57834",
    "CHEBI:58704": "CHEBI:57920",
}


class SpermineOxidase(ReactionClass):
    """spermine oxidase activity.

    Catalysis of the reaction: H2O + O2 + spermine = 3-aminopropanal + H2O2 + spermidine. Weak activity with N(1)-acetylspermine. The Arabidopsis thaliana enzyme converts norspermine to norspermidine.
    """

    GO_ID = "GO:0052901"
    EC_NUMBER_PREFIX = "1.5.3.16"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}

        if CHEBI_O2 not in left_chebis or CHEBI_H2O not in left_chebis:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and water substrates")
        if CHEBI_H2O2 not in right_chebis or CHEBI_AMINOPROPANAL not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Requires hydrogen peroxide and 3-aminopropanal products",
            )

        substrate = next((chebi for chebi in POLYAMINE_PRODUCT_PAIRS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(
                is_member=False,
                explanation="No spermine or norspermine substrate detected",
            )

        product = POLYAMINE_PRODUCT_PAIRS[substrate]
        if product not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing the expected shorter polyamine product",
            )

        substantive_left = left_chebis - {CHEBI_O2, CHEBI_H2O, substrate}
        substantive_right = right_chebis - {CHEBI_H2O2, CHEBI_AMINOPROPANAL, product}
        if substantive_left or substantive_right:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside spermine oxidase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Spermine oxidase: oxygen-dependent oxidative cleavage of spermine or norspermine to 3-aminopropanal and a shorter polyamine",
        )
