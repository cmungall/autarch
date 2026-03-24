"""deoxynucleoside kinase activity.

Catalysis of ATP-dependent phosphorylation of deoxynucleosides to the
corresponding 5'-monophosphates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class DeoxynucleosideKinase(ReactionClass):
    """deoxynucleoside kinase activity.

    Catalysis of ATP-dependent phosphorylation of deoxynucleosides to the
    corresponding 5'-monophosphates.
    """

    GO_ID = "GO:0019136"
    EC_NUMBER_PREFIX = "2.7.1.145"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:18274": "CHEBI:65317",
        "CHEBI:17748": "CHEBI:63528",
        "CHEBI:17172": "CHEBI:57673",
        "CHEBI:17256": "CHEBI:58245",
        "CHEBI:27596": "CHEBI:58777",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_ATP not in left_ids or CHEBI_ADP not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires ATP donor and ADP product",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_ATP, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_ADP, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one deoxynucleoside substrate and one monophosphate product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(
                is_member=False,
                explanation="Substrate/product branch does not match a supported deoxynucleoside kinase pair",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Deoxynucleoside kinase: ATP-dependent phosphorylation of a deoxynucleoside",
        )
