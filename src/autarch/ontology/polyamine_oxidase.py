"""polyamine oxidase activity.

Catalysis of oxidation of polyamines with dioxygen, producing hydrogen peroxide
and a shortened polyamine or aminoaldehyde product set.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_H2O2 = "CHEBI:16240"


class PolyamineOxidase(ReactionClass):
    """polyamine oxidase activity.

    Catalysis of oxidation of polyamines with dioxygen, producing hydrogen peroxide
    and a shortened polyamine or aminoaldehyde product set.
    """

    GO_ID = "GO:0046592"
    EC_NUMBER_PREFIX = "1.5.3.17"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:58101": {"CHEBI:30322", "CHEBI:57834"},
        "CHEBI:45725": {"CHEBI:133427", "CHEBI:57834"},
        "CHEBI:57834": {"CHEBI:133427", "CHEBI:326268"},
        "CHEBI:58324": {"CHEBI:30322", "CHEBI:326268"},
        "CHEBI:58535": {"CHEBI:7386", "CHEBI:57484"},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_O2 not in left_ids or CHEBI_H2O not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires dioxygen and water substrates",
            )
        if CHEBI_H2O2 not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires hydrogen peroxide product",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_O2, CHEBI_H2O}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_H2O2]
        if len(left_core) != 1 or len(right_core) != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one polyamine substrate and two substantive products",
            )

        substrate = left_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCTS:
            return ClassificationResult(
                is_member=False,
                explanation="No supported polyamine substrate detected",
            )
        if set(right_core) != self.SUBSTRATE_TO_PRODUCTS[substrate]:
            return ClassificationResult(
                is_member=False,
                explanation="Products do not match a supported polyamine oxidase branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Polyamine oxidase: oxidative cleavage of a polyamine with hydrogen peroxide formation",
        )
