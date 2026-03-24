"""nucleotide diphosphatase.

Catalysis of the reaction: a nucleoside triphosphate + H2O = the corresponding
nucleoside monophosphate + diphosphate + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_DIPHOSPHATE = "CHEBI:33019"


class NucleotideDiphosphatase(ReactionClass):
    """nucleotide diphosphatase.

    Catalysis of the reaction: a nucleoside triphosphate + H2O = the corresponding
    nucleoside monophosphate + diphosphate + H(+).
    """

    EC_NUMBER_PREFIX = "3.6.1.9"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:61555": "CHEBI:246422",
        "CHEBI:30616": "CHEBI:456215",
        "CHEBI:61481": "CHEBI:57566",
        "CHEBI:61557": "CHEBI:58043",
        "CHEBI:37563": "CHEBI:60377",
        "CHEBI:61404": "CHEBI:58245",
        "CHEBI:61382": "CHEBI:61194",
        "CHEBI:61429": "CHEBI:57673",
        "CHEBI:37568": "CHEBI:63528",
        "CHEBI:61314": "CHEBI:57464",
        "CHEBI:37565": "CHEBI:58115",
        "CHEBI:46398": "CHEBI:57865",
        "CHEBI:61402": "CHEBI:58053",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_H2O not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires water as reactant")
        if CHEBI_DIPHOSPHATE not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires diphosphate product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_DIPHOSPHATE, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one nucleotide substrate and one nucleotide monophosphate product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported nucleotide diphosphatase pair")

        return ClassificationResult(is_member=True, explanation="Nucleotide diphosphatase: hydrolysis of a nucleoside triphosphate to the corresponding monophosphate and diphosphate")
