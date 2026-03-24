"""aromatic-amino-acid transaminase.

Catalysis of transamination between 2-oxoglutarate and supported amino-acid
substrates to produce L-glutamate and the corresponding oxo-acid product.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_GLUTAMATE = "CHEBI:29985"


class AromaticAminoAcidTransaminase(ReactionClass):
    """aromatic-amino-acid transaminase.

    Catalysis of transamination between 2-oxoglutarate and supported amino-acid
    substrates to produce L-glutamate and the corresponding oxo-acid product.
    """

    EC_NUMBER_PREFIX = "2.6.1.57"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:58315": "CHEBI:36242",
        "CHEBI:58095": "CHEBI:18005",
        "CHEBI:57844": "CHEBI:133493",
        "CHEBI:84824": "CHEBI:73309",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_2_OXOGLUTARATE not in left_ids or CHEBI_GLUTAMATE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires 2-oxoglutarate substrate and glutamate product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_2_OXOGLUTARATE]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_GLUTAMATE]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one amino-acid substrate and one oxo-acid product")

        if self.SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Substrate/product branch does not match a supported aromatic-amino-acid transaminase pair")

        return ClassificationResult(is_member=True, explanation="Aromatic-amino-acid transaminase: 2-oxoglutarate-dependent transamination of a supported substrate")
