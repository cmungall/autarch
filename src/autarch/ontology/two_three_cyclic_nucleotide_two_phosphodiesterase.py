"""2',3'-cyclic-nucleotide 2'-phosphodiesterase.

Catalysis of the reaction: a nucleoside 2',3'-cyclic phosphate + H2O = the
corresponding nucleoside 3'-phosphate + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class TwoThreeCyclicNucleotideTwoPhosphodiesterase(ReactionClass):
    """2',3'-cyclic-nucleotide 2'-phosphodiesterase.

    Catalysis of the reaction: a nucleoside 2',3'-cyclic phosphate + H2O = the
    corresponding nucleoside 3'-phosphate + H(+).
    """

    EC_NUMBER_PREFIX = "3.1.4.16"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:66954": "CHEBI:66949",
        "CHEBI:60837": "CHEBI:60732",
        "CHEBI:60873": "CHEBI:60784",
        "CHEBI:60877": "CHEBI:60875",
        "CHEBI:60879": "CHEBI:60880",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires water as reactant")
        if CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_H_PLUS]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one cyclic nucleotide substrate and one 3'-phosphate product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported 2',3'-cyclic-nucleotide 2'-phosphodiesterase pair")

        return ClassificationResult(is_member=True, explanation="2',3'-cyclic-nucleotide 2'-phosphodiesterase: hydrolysis of a cyclic nucleotide to the corresponding 3'-phosphate")
