"""nucleoside diphosphate kinase activity.

Catalysis of the reaction: ATP + nucleoside diphosphate = ADP + nucleoside
triphosphate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class NucleosideDiphosphateKinase(ReactionClass):
    """nucleoside diphosphate kinase activity.

    Catalysis of the reaction: ATP + nucleoside diphosphate = ADP + nucleoside
    triphosphate.
    """

    GO_ID = "GO:0004550"
    EC_NUMBER_PREFIX = "2.7.4.6"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:57930": "CHEBI:61557",
        "CHEBI:58223": "CHEBI:46398",
        "CHEBI:58069": "CHEBI:37563",
        "CHEBI:57667": "CHEBI:61404",
        "CHEBI:58593": "CHEBI:61481",
        "CHEBI:58369": "CHEBI:37568",
        "CHEBI:58189": "CHEBI:37565",
        "CHEBI:58595": "CHEBI:61429",
        "CHEBI:60471": "CHEBI:61555",
        "CHEBI:58280": "CHEBI:61402",
        "CHEBI:62286": "CHEBI:61382",
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
        if CHEBI_ATP not in left_ids or CHEBI_ADP not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires ATP donor and ADP product",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_ATP]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_ADP]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one nucleoside diphosphate substrate and one nucleoside triphosphate product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCT:
            return ClassificationResult(
                is_member=False,
                explanation="No supported nucleoside diphosphate substrate detected",
            )
        if self.SUBSTRATE_TO_PRODUCT[substrate] != product:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not match the expected nucleoside triphosphate branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Nucleoside diphosphate kinase: ATP-dependent phosphorylation of a nucleoside diphosphate",
        )
