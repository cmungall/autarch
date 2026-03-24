"""short-chain 2-methyl fatty acyl-CoA dehydrogenase activity.

Catalysis of electron-transfer-flavoprotein-linked oxidation of short-chain
2-methyl fatty acyl-CoA substrates to the corresponding enoyl-CoA products.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_OXIDIZED_ETF = "CHEBI:57692"
CHEBI_REDUCED_ETF = "CHEBI:58307"


class ShortChain2MethylFattyAcylCoADehydrogenaseActivity(ReactionClass):
    """short-chain 2-methyl fatty acyl-CoA dehydrogenase activity."""

    GO_ID = "GO:0003853"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:57336": "CHEBI:57337",
        "CHEBI:57338": "CHEBI:62500",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            substrate_to_product=self.SUBSTRATE_TO_PRODUCT,
            explanation="short-chain 2-methyl fatty acyl-CoA dehydrogenase activity: ETF-linked oxidation of a short-chain 2-methyl acyl-CoA",
        )
        if forward.is_member:
            return forward
        reverse = self._check_direction(
            left_ids=[p.chebi_id for p in reaction.right_participants if p.chebi_id],
            right_ids=[p.chebi_id for p in reaction.left_participants if p.chebi_id],
            substrate_to_product={v: k for k, v in self.SUBSTRATE_TO_PRODUCT.items()},
            explanation="short-chain 2-methyl fatty acyl-CoA dehydrogenase activity: ETF-linked oxidation of a short-chain 2-methyl acyl-CoA (reverse reaction orientation)",
        )
        if reverse.is_member:
            return reverse
        return forward

    def _check_direction(
        self,
        left_ids: list[str],
        right_ids: list[str],
        substrate_to_product: dict[str, str],
        explanation: str,
    ) -> ClassificationResult:
        if CHEBI_OXIDIZED_ETF not in left_ids or CHEBI_REDUCED_ETF not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires oxidized/reduced electron-transfer flavoprotein pair")
        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_OXIDIZED_ETF, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {CHEBI_REDUCED_ETF, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one acyl-CoA substrate and one enoyl-CoA product")
        if substrate_to_product.get(left_core[0]) != right_core[0]:
            return ClassificationResult(is_member=False, explanation="Acyl-CoA branch does not match a supported short-chain 2-methyl substrate pair")
        return ClassificationResult(is_member=True, explanation=explanation)
