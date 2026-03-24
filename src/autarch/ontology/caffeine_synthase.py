"""caffeine synthase activity.

Catalysis of the reaction: 7-methylxanthine/theobromine/1,7-dimethylxanthine + S-adenosyl-L-methionine = methylated xanthine + S-adenosyl-L-homocysteine + H(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

METHYLATION_PAIRS = {
    "CHEBI:48991": "CHEBI:28946",
    "CHEBI:28946": "CHEBI:27732",
    "CHEBI:25858": "CHEBI:27732",
}


class CaffeineSynthase(Methyltransferase):
    """caffeine synthase activity.

    Catalysis of the reaction: 7-methylxanthine/theobromine/1,7-dimethylxanthine + S-adenosyl-L-methionine = methylated xanthine + S-adenosyl-L-homocysteine + H(+).
    """

    GO_ID = "GO:0102741"
    EC_NUMBER_PREFIX = "2.1.1.160"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one xanthine substrate and one methylated xanthine product")

        substrate = left_core[0].chebi_id
        product = right_core[0].chebi_id
        if substrate not in METHYLATION_PAIRS or METHYLATION_PAIRS[substrate] != product:
            return ClassificationResult(is_member=False, explanation="Requires one of the known caffeine-synthase xanthine methylation steps")

        return ClassificationResult(
            is_member=True,
            explanation="Caffeine synthase: SAM-dependent methylation within the methylxanthine pathway to caffeine",
        )
