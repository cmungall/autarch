"""2'-deamino-2'-hydroxyneamine 1-alpha-D-kanosaminyltransferase.

Catalysis of transfer of alpha-D-kanosamine from UDP-alpha-D-kanosamine to
supported aminoglycoside acceptor substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_UDP = "CHEBI:58223"
CHEBI_UDP_ALPHA_D_KANOSAMINE = "CHEBI:71964"


class TwoDeaminoTwoHydroxyneamine1AlphaDKanosaminyltransferase(ReactionClass):
    """2'-deamino-2'-hydroxyneamine 1-alpha-D-kanosaminyltransferase.

    Catalysis of transfer of alpha-D-kanosamine from UDP-alpha-D-kanosamine to
    supported aminoglycoside acceptor substrates.
    """

    EC_NUMBER_PREFIX = "2.4.1.301"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:65076": "CHEBI:58549",
        "CHEBI:65015": "CHEBI:72755",
        "CHEBI:65071": "CHEBI:72756",
        "CHEBI:67213": "CHEBI:58214",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_UDP_ALPHA_D_KANOSAMINE not in left_ids or CHEBI_UDP not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires UDP-alpha-D-kanosamine donor and UDP coproduct")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_UDP_ALPHA_D_KANOSAMINE, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_UDP, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one aminoglycoside acceptor and one kanosaminylated product")

        if self.SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Aminoglycoside branch does not match a supported kanosaminyltransferase pair")

        return ClassificationResult(is_member=True, explanation="2'-deamino-2'-hydroxyneamine 1-alpha-D-kanosaminyltransferase: transfer of alpha-D-kanosamine from UDP-alpha-D-kanosamine")
