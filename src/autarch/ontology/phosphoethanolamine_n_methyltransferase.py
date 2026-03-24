"""phosphoethanolamine N-methyltransferase activity.

Catalysis of the reaction: phosphoethanolamine/N-methylethanolamine phosphate/N,N-dimethylethanolamine phosphate + S-adenosyl-L-methionine = methylated phosphobase + S-adenosyl-L-homocysteine + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

METHYLATION_PAIRS = {
    "CHEBI:58190": "CHEBI:57781",
    "CHEBI:57781": "CHEBI:58641",
    "CHEBI:58641": "CHEBI:295975",
}


class PhosphoethanolamineNMethyltransferase(Methyltransferase):
    """phosphoethanolamine N-methyltransferase activity.

    Catalysis of the reaction: phosphoethanolamine/N-methylethanolamine phosphate/N,N-dimethylethanolamine phosphate + S-adenosyl-L-methionine = methylated phosphobase + S-adenosyl-L-homocysteine + H(+).
    """

    GO_ID = "GO:0000234"
    EC_NUMBER_PREFIX = "2.1.1.103"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one phosphobase substrate and one methylated phosphobase product")

        substrate = left_core[0].chebi_id
        product = right_core[0].chebi_id
        if substrate not in METHYLATION_PAIRS or METHYLATION_PAIRS[substrate] != product:
            return ClassificationResult(is_member=False, explanation="Requires one of the phosphoethanolamine N-methyltransferase methylation steps")

        return ClassificationResult(
            is_member=True,
            explanation="Phosphoethanolamine N-methyltransferase: SAM-dependent stepwise N-methylation of phosphoethanolamine derivatives",
        )
