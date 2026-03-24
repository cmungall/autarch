"""D-aminoacyl-tRNA deacylase activity.

Catalysis of the reaction: a D-aminoacyl-tRNA + H2O = a D-alpha-amino acid + a
 tRNA + H(+). Removal of a D-amino acid from a charged tRNA.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_TRNA = "CHEBI:78442"


class DAminoacylTRNADeacylase(ReactionClass):
    """D-aminoacyl-tRNA deacylase activity.

    Catalysis of the reaction: a D-aminoacyl-tRNA + H2O = a D-alpha-amino acid + a
    tRNA + H(+). Removal of a D-amino acid from a charged tRNA.
    """

    GO_ID = "GO:0051499"
    EC_NUMBER_PREFIX = "3.1.1.96"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    CHARGED_TO_FREE = {
        "CHEBI:79333": "CHEBI:59871",
        "CHEBI:78723": "CHEBI:58570",
        "CHEBI:188448": "CHEBI:57981",
        "CHEBI:188450": "CHEBI:57719",
        "CHEBI:188449": "CHEBI:29990",
        "CHEBI:78522": "CHEBI:57305",
    }
    FREE_TO_CHARGED = {value: key for key, value in CHARGED_TO_FREE.items()}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_forward(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        )
        if forward.is_member:
            return forward

        reverse = self._check_reverse(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
        )
        if reverse.is_member:
            return reverse

        return forward

    def _check_forward(self, left_ids: list[str], right_ids: list[str]) -> ClassificationResult:
        if CHEBI_H2O not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires water for D-aminoacyl-tRNA hydrolysis",
            )

        charged = [chebi_id for chebi_id in left_ids if chebi_id != CHEBI_H2O]
        products = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_H_PLUS]
        if len(charged) != 1 or len(products) != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one D-aminoacyl-tRNA substrate and two core products",
            )

        substrate = charged[0]
        if substrate not in self.CHARGED_TO_FREE:
            return ClassificationResult(
                is_member=False,
                explanation="No supported D-aminoacyl-tRNA substrate detected",
            )
        if set(products) != {self.CHARGED_TO_FREE[substrate], CHEBI_TRNA}:
            return ClassificationResult(
                is_member=False,
                explanation="Products do not match the expected D-amino acid and tRNA pair",
            )

        return ClassificationResult(
            is_member=True,
            explanation="D-aminoacyl-tRNA deacylase: hydrolysis of a D-aminoacyl-tRNA to free D-amino acid and tRNA",
        )

    def _check_reverse(self, left_ids: list[str], right_ids: list[str]) -> ClassificationResult:
        if CHEBI_H2O not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Reverse orientation must still include water as the byproduct",
            )

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {CHEBI_TRNA, CHEBI_H_PLUS}]
        if CHEBI_TRNA not in left_ids or len(left_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected free D-amino acid plus tRNA on the left in reverse orientation",
            )

        right_core = [chebi_id for chebi_id in right_ids if chebi_id != CHEBI_H2O]
        if len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one D-aminoacyl-tRNA product in reverse orientation",
            )

        substrate = left_core[0]
        product = right_core[0]
        if self.FREE_TO_CHARGED.get(substrate) != product:
            return ClassificationResult(
                is_member=False,
                explanation="Reverse orientation does not match a supported D-aminoacyl-tRNA branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="D-aminoacyl-tRNA deacylase: reverse reaction orientation for a supported D-aminoacyl-tRNA branch",
        )
