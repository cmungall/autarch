"""gibberellin 20-oxidase activity.

Catalysis of the reaction: 2 2-oxoglutarate + gibberellin A12 (GA12) + H+ + 3
O2 = 3 CO2 + gibberellin A9 (GA9) + 2 H2O + 2 succinate. This reaction results
in the oxidation of C-20 gibberellins to form the corresponding C-19 lactones,
via a three-step oxidation at C-20 of the GA skeleton. Also converts GA53 to
GA20. GA25 is also formed as a minor product.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_GA12 = "CHEBI:58627"
CHEBI_GA15 = "CHEBI:143956"
CHEBI_GA24 = "CHEBI:143957"
CHEBI_GA25 = "CHEBI:143959"
CHEBI_GA9 = "CHEBI:73255"
CHEBI_GA53 = "CHEBI:143954"
CHEBI_GA44 = "CHEBI:143955"
CHEBI_GA20 = "CHEBI:58526"


class Gibberellin20Oxidase(ReactionClass):
    """gibberellin 20-oxidase activity.

    Catalysis of the reaction: 2 2-oxoglutarate + gibberellin A12 (GA12) +
    H+ + 3 O2 = 3 CO2 + gibberellin A9 (GA9) + 2 H2O + 2 succinate. This
    reaction results in the oxidation of C-20 gibberellins to form the
    corresponding C-19 lactones, via a three-step oxidation at C-20 of the GA
    skeleton. Also converts GA53 to GA20. GA25 is also formed as a minor
    product.
    """

    GO_ID = "GO:0045544"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        CHEBI_GA12: {CHEBI_GA25, CHEBI_GA15, CHEBI_GA9},
        CHEBI_GA15: {CHEBI_GA24},
        CHEBI_GA53: {CHEBI_GA44, CHEBI_GA20},
    }
    LEFT_SPECTATORS = {CHEBI_2_OXOGLUTARATE, CHEBI_O2, CHEBI_H_PLUS}
    RIGHT_SPECTATORS = {CHEBI_SUCCINATE, CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_2_OXOGLUTARATE not in left_ids or CHEBI_O2 not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires 2-oxoglutarate and dioxygen cosubstrates",
            )
        if CHEBI_SUCCINATE not in right_ids or CHEBI_CO2 not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires succinate and carbon dioxide coproducts",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in self.LEFT_SPECTATORS]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in self.RIGHT_SPECTATORS]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one gibberellin substrate and one oxidized gibberellin product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCTS:
            return ClassificationResult(
                is_member=False,
                explanation="No supported gibberellin 20-oxidase substrate detected",
            )
        if product not in self.SUBSTRATE_TO_PRODUCTS[substrate]:
            return ClassificationResult(
                is_member=False,
                explanation="Product is not a recognized gibberellin 20-oxidase oxidation branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Gibberellin 20-oxidase: 2-oxoglutarate-dependent oxidation of a C-20 gibberellin substrate",
        )
