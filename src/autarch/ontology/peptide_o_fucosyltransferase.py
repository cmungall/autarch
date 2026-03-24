"""peptide-O-fucosyltransferase activity.

Catalysis of transfer of alpha-L-fucose from GDP-beta-L-fucose to serine or
threonine residues in peptide or protein substrates.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_GDP = "CHEBI:58189"
CHEBI_GDP_BETA_L_FUCOSE = "CHEBI:57273"


class PeptideOFucosyltransferase(ReactionClass):
    """peptide-O-fucosyltransferase activity.

    Catalysis of transfer of alpha-L-fucose from GDP-beta-L-fucose to serine or
    threonine residues in peptide or protein substrates.
    """

    GO_ID = "GO:0046922"
    EC_NUMBER_PREFIX = "2.4.1.221"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:29999": "CHEBI:189632",
        "CHEBI:30013": "CHEBI:189631",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            reverse=True,
        )
        if reverse.is_member:
            return reverse

        return forward

    def _check_direction(self, left_ids: list[str], right_ids: list[str], reverse: bool = False) -> ClassificationResult:
        donor = CHEBI_GDP if reverse else CHEBI_GDP_BETA_L_FUCOSE
        coproduct = CHEBI_GDP_BETA_L_FUCOSE if reverse else CHEBI_GDP
        if donor not in left_ids or coproduct not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires GDP-beta-L-fucose/GDP transfer pair",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {donor, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {coproduct, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one peptide substrate and one fucosylated peptide product",
            )

        mapping = {v: k for k, v in self.SUBSTRATE_TO_PRODUCT.items()} if reverse else self.SUBSTRATE_TO_PRODUCT
        if mapping.get(left_core[0]) != right_core[0]:
            return ClassificationResult(
                is_member=False,
                explanation="Peptide branch does not match a supported O-fucosylation pair",
            )

        detail = "reverse reaction orientation" if reverse else "GDP-fucose-dependent O-fucosylation of a peptide residue"
        return ClassificationResult(is_member=True, explanation=f"Peptide-O-fucosyltransferase: {detail}")
