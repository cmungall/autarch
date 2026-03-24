"""hexokinase.

Catalysis of the reaction: a D-hexose + ATP = the corresponding D-hexose
6-phosphate + ADP + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_D_GLUCOSAMINE = "CHEBI:58723"
CHEBI_D_GLUCOSAMINE_6_PHOSPHATE = "CHEBI:58725"
CHEBI_D_MANNOSE = "CHEBI:4208"
CHEBI_D_MANNOSE_6_PHOSPHATE = "CHEBI:58735"
CHEBI_D_FRUCTOSE = "CHEBI:37721"
CHEBI_D_FRUCTOSE_6_PHOSPHATE = "CHEBI:61527"
CHEBI_D_GLUCOSE = "CHEBI:4167"
CHEBI_D_GLUCOSE_6_PHOSPHATE = "CHEBI:61548"
CHEBI_D_HEXOSE = "CHEBI:4194"
CHEBI_D_HEXOSE_6_PHOSPHATE = "CHEBI:229467"


class Hexokinase(ReactionClass):
    """hexokinase.

    Catalysis of the reaction: a D-hexose + ATP = the corresponding D-hexose
    6-phosphate + ADP + H(+).
    """

    EC_NUMBER_PREFIX = "2.7.1.1"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        CHEBI_D_GLUCOSAMINE: CHEBI_D_GLUCOSAMINE_6_PHOSPHATE,
        CHEBI_D_MANNOSE: CHEBI_D_MANNOSE_6_PHOSPHATE,
        CHEBI_D_FRUCTOSE: CHEBI_D_FRUCTOSE_6_PHOSPHATE,
        CHEBI_D_GLUCOSE: CHEBI_D_GLUCOSE_6_PHOSPHATE,
        CHEBI_D_HEXOSE: CHEBI_D_HEXOSE_6_PHOSPHATE,
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
                explanation="Expected one hexose substrate and one hexose 6-phosphate product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if substrate not in self.SUBSTRATE_TO_PRODUCT:
            return ClassificationResult(
                is_member=False,
                explanation="No supported D-hexose substrate detected",
            )
        if self.SUBSTRATE_TO_PRODUCT[substrate] != product:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not match the expected hexose 6-phosphate branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Hexokinase: ATP-dependent phosphorylation of a D-hexose at the 6-position",
        )
