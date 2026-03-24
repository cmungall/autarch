"""N(1)-acetylpolyamine oxidase (3-acetamidopropanal-forming) activity.

Catalysis of the reaction: H2O + N(1)-acetylspermine + O2 = 3-acetamidopropanal + H2O2 + spermidine. Also converts N(1)-acetylspermidine to putrescine.
"""

from abc import abstractmethod

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H2O2, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_N1_ACETYLSPERMINE = "CHEBI:58101"
CHEBI_N1_ACETYLSPERMIDINE = "CHEBI:58324"
CHEBI_DIACETYLSPERMINE = "CHEBI:58550"
CHEBI_3_ACETAMIDOPROPANAL = "CHEBI:30322"
CHEBI_SPERMIDINE = "CHEBI:57834"
CHEBI_PUTRESCINE = "CHEBI:326268"

PRODUCT_PAIRS = {
    CHEBI_N1_ACETYLSPERMINE: CHEBI_SPERMIDINE,
    CHEBI_N1_ACETYLSPERMIDINE: CHEBI_PUTRESCINE,
    CHEBI_DIACETYLSPERMINE: CHEBI_N1_ACETYLSPERMIDINE,
}


class N1AcetylpolyamineOxidase3AcetamidopropanalForming(ReactionClass):
    """Legacy implementation for N(1)-acetylpolyamine oxidase.

    Catalysis of the reaction: H2O + N(1)-acetylspermine + O2 = 3-acetamidopropanal + H2O2 + spermidine. Also converts N(1)-acetylspermidine to putrescine.
    """

    EC_NUMBER_PREFIX = "1.5.3.13"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    @abstractmethod
    def _implementation_only(self) -> None:
        """Mark this legacy implementation class as abstract."""

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if CHEBI_O2 not in left_chebis or CHEBI_H2O not in left_chebis:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and water substrates")
        if CHEBI_H2O2 not in right_chebis or CHEBI_3_ACETAMIDOPROPANAL not in right_chebis:
            return ClassificationResult(is_member=False, explanation="Missing hydrogen peroxide or 3-acetamidopropanal product")

        substrate = next((chebi for chebi in PRODUCT_PAIRS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(is_member=False, explanation="No N(1)-acetylpolyamine substrate detected")

        expected_polyamine_product = PRODUCT_PAIRS[substrate]
        if expected_polyamine_product not in right_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing the expected deacetylated polyamine product",
            )

        substantive_left = left_chebis - {CHEBI_O2, CHEBI_H2O}
        substantive_right = right_chebis - {CHEBI_H2O2, CHEBI_3_ACETAMIDOPROPANAL}
        if substantive_left != {substrate} or substantive_right != {expected_polyamine_product}:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside N(1)-acetylpolyamine oxidase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="N(1)-acetylpolyamine oxidase: oxygen-dependent oxidative cleavage of an acetylpolyamine to 3-acetamidopropanal and a shorter polyamine",
        )
