"""purine-nucleoside phosphorylase activity.

Catalysis of the reaction: purine nucleoside + phosphate = purine + alpha-D-ribose
1-phosphate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PHOSPHATE = "CHEBI:16838"
CHEBI_RIBOSE_1_PHOSPHATE = "CHEBI:57720"
CHEBI_DEOXY_RIBOSE_1_PHOSPHATE = "CHEBI:57259"


class PurineNucleosidePhosphorylase(ReactionClass):
    """purine-nucleoside phosphorylase activity.

    Catalysis of the reaction: purine nucleoside + phosphate = purine + alpha-D-ribose
    1-phosphate.
    """

    GO_ID = "GO:0004731"
    EC_NUMBER_PREFIX = "2.4.2.1"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCTS = {
        "CHEBI:16750": {"CHEBI:16235", CHEBI_RIBOSE_1_PHOSPHATE},
        "CHEBI:18107": {"CHEBI:17712", CHEBI_RIBOSE_1_PHOSPHATE},
        "CHEBI:16335": {"CHEBI:16708", CHEBI_RIBOSE_1_PHOSPHATE},
        "CHEBI:17596": {"CHEBI:17368", CHEBI_RIBOSE_1_PHOSPHATE},
        "CHEBI:17172": {"CHEBI:16235", CHEBI_DEOXY_RIBOSE_1_PHOSPHATE},
        "CHEBI:17256": {"CHEBI:16708", CHEBI_DEOXY_RIBOSE_1_PHOSPHATE},
        "CHEBI:28997": {"CHEBI:17368", CHEBI_DEOXY_RIBOSE_1_PHOSPHATE},
        "CHEBI:142355": {"CHEBI:26386", CHEBI_RIBOSE_1_PHOSPHATE},
        "CHEBI:142361": {"CHEBI:26386", CHEBI_DEOXY_RIBOSE_1_PHOSPHATE},
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
        )
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_direction(self, left_ids: list[str], right_ids: list[str]) -> ClassificationResult:
        if CHEBI_PHOSPHATE not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires phosphate as a cosubstrate",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_PHOSPHATE]
        if len(left_core) != 1 or len(right_ids) != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one purine nucleoside substrate and two products",
            )

        substrate = left_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCTS:
            return ClassificationResult(
                is_member=False,
                explanation="No supported purine nucleoside substrate detected",
            )

        if set(right_ids) != self.SUBSTRATE_TO_PRODUCTS[substrate]:
            return ClassificationResult(
                is_member=False,
                explanation="Products do not match the expected purine base and ribose 1-phosphate branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Purine-nucleoside phosphorylase: phosphorolysis of a purine nucleoside to a purine base and ribose 1-phosphate",
        )
