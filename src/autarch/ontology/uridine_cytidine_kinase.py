"""uridine/cytidine kinase.

Catalysis of nucleoside-triphosphate-dependent phosphorylation of uridine or
cytidine to the corresponding 5'-monophosphates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_GDP, CHEBI_GTP, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class UridineCytidineKinase(ReactionClass):
    """uridine/cytidine kinase.

    Catalysis of nucleoside-triphosphate-dependent phosphorylation of uridine or
    cytidine to the corresponding 5'-monophosphates.
    """

    EC_NUMBER_PREFIX = "2.7.1.48"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:16704": "CHEBI:57865",
        "CHEBI:17562": "CHEBI:60377",
    }
    DONOR_TO_COPRODUCT = {
        CHEBI_ATP: CHEBI_ADP,
        CHEBI_GTP: CHEBI_GDP,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        matched_pair = next(((donor, coproduct) for donor, coproduct in self.DONOR_TO_COPRODUCT.items() if donor in left_ids and coproduct in right_ids), None)
        if matched_pair is None:
            return ClassificationResult(is_member=False, explanation="Requires ATP or GTP donor with matching diphosphate coproduct")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {matched_pair[0], CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {matched_pair[1], CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one nucleoside substrate and one monophosphate product")

        if self.SUBSTRATE_TO_PRODUCT.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported uridine/cytidine kinase pair")

        return ClassificationResult(is_member=True, explanation="Uridine/cytidine kinase: nucleoside-triphosphate-dependent phosphorylation of uridine or cytidine")
