"""(R)-amidase.

Catalysis of hydrolysis of supported chiral amides to the corresponding acids or
amino-carboxylates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_AMMONIUM = "CHEBI:28938"


class RAmidase(ReactionClass):
    """(R)-amidase.

    Catalysis of hydrolysis of supported chiral amides to the corresponding acids or
    amino-carboxylates.
    """

    EC_NUMBER_PREFIX = "3.5.1.100"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:58916": {"CHEBI:58917", CHEBI_AMMONIUM},
        "CHEBI:58918": {"CHEBI:57966", CHEBI_AMMONIUM},
        "CHEBI:60118": {"CHEBI:60120", CHEBI_AMMONIUM},
        "CHEBI:60254": {"CHEBI:224366", "CHEBI:58917"},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires water substrate")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        if len(left_core) != 1 or len(right_ids) != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one amide substrate and two products",
            )

        substrate = left_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCTS:
            return ClassificationResult(is_member=False, explanation="No supported (R)-amidase substrate detected")
        if set(right_ids) != self.SUBSTRATE_TO_PRODUCTS[substrate]:
            return ClassificationResult(is_member=False, explanation="Products do not match a supported (R)-amidase hydrolysis branch")

        return ClassificationResult(is_member=True, explanation="(R)-amidase: hydrolysis of a supported chiral amide substrate")
