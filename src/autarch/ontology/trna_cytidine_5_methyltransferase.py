"""tRNA (cytidine-5-)-methyltransferase activity.

Catalysis of SAM-dependent methylation of cytidine residues in tRNA or tRNA
precursor molecules to 5-methylcytidine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_CYTIDINE_IN_TRNA = "CHEBI:82748"
CHEBI_METHYLCYTIDINE_IN_TRNA = "CHEBI:74483"


class TRNACytidine5Methyltransferase(ReactionClass):
    """tRNA (cytidine-5-)-methyltransferase activity.

    Catalysis of SAM-dependent methylation of cytidine residues in tRNA or tRNA
    precursor molecules to 5-methylcytidine.
    """

    GO_ID = "GO:0016428"
    EC_NUMBER_PREFIX = "2.1.1.202"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_sam_exact_pair(
            reaction,
            substrate=CHEBI_CYTIDINE_IN_TRNA,
            product=CHEBI_METHYLCYTIDINE_IN_TRNA,
            explanation="tRNA (cytidine-5-)-methyltransferase: SAM-dependent methylation of cytidine in tRNA",
        )


def _check_reversible_sam_exact_pair(
    reaction: Reaction,
    substrate: str,
    product: str,
    explanation: str,
) -> ClassificationResult:
    forward = _check_sam_direction(
        left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
        right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        substrate=substrate,
        product=product,
        donor=CHEBI_SAM,
        coproduct=CHEBI_SAH,
        explanation=explanation,
    )
    if forward.is_member:
        return forward

    reverse = _check_sam_direction(
        left_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        right_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
        substrate=product,
        product=substrate,
        donor=CHEBI_SAH,
        coproduct=CHEBI_SAM,
        explanation=f"{explanation} (reverse reaction orientation)",
    )
    if reverse.is_member:
        return reverse

    return forward


def _check_sam_direction(
    left_ids: list[str],
    right_ids: list[str],
    substrate: str,
    product: str,
    donor: str,
    coproduct: str,
    explanation: str,
) -> ClassificationResult:
    if donor not in left_ids or coproduct not in right_ids:
        return ClassificationResult(is_member=False, explanation="Requires SAM/SAH cofactor pair")

    left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {donor, CHEBI_H_PLUS}]
    right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {coproduct, CHEBI_H_PLUS}]
    if len(left_core) != 1 or len(right_core) != 1:
        return ClassificationResult(is_member=False, explanation="Expected one RNA substrate and one methylated RNA product")

    if left_core[0] != substrate or right_core[0] != product:
        return ClassificationResult(is_member=False, explanation="RNA branch does not match a supported methylation pair")

    return ClassificationResult(is_member=True, explanation=explanation)
