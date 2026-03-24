"""tocopherol C-methyltransferase activity.

Catalysis of the reaction: gamma-tocopherol + S-adenosyl-L-methionine = (+)-alpha-tocopherol + H+ + S-adenosyl-L-homocysteine. This reaction can also use delta-tocopherol, gamma-tocotrienol and delta-tocotrienol as substrates.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_GAMMA_TOCOPHEROL = "CHEBI:18185"
CHEBI_ALPHA_TOCOPHEROL = "CHEBI:18145"


class TocopherolCMethyltransferase(Methyltransferase):
    """tocopherol C-methyltransferase activity.

    Catalysis of the reaction: gamma-tocopherol + S-adenosyl-L-methionine = (+)-alpha-tocopherol + H+ + S-adenosyl-L-homocysteine. This reaction can also use delta-tocopherol, gamma-tocotrienol and delta-tocotrienol as substrates.
    """

    GO_ID = "GO:0050342"
    EC_NUMBER_PREFIX = "2.1.1.95"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [p for p in reaction.left_participants if p.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}]
        right_core = [p for p in reaction.right_participants if p.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one tocopherol substrate/product pair")

        substrate = left_core[0]
        product = right_core[0]
        if substrate.chebi_id == CHEBI_GAMMA_TOCOPHEROL and product.chebi_id == CHEBI_ALPHA_TOCOPHEROL:
            return ClassificationResult(
                is_member=True,
                explanation="Tocopherol C-methyltransferase: SAM-dependent methylation of a tocopherol scaffold to alpha-tocopherol",
            )

        if not (self._is_tocopherol_like(substrate) and self._is_tocopherol_like(product)):
            return ClassificationResult(is_member=False, explanation="No tocopherol-like substrate/product pair detected")
        if self._carbon_count(product) != self._carbon_count(substrate) + 1:
            return ClassificationResult(is_member=False, explanation="Product is not a single-carbon methyl homolog of the tocopherol substrate")
        if self._oxygen_count(product) != self._oxygen_count(substrate):
            return ClassificationResult(is_member=False, explanation="Tocopherol oxygen count is not preserved across methyl transfer")

        return ClassificationResult(
            is_member=True,
            explanation="Tocopherol C-methyltransferase: SAM-dependent carbon methylation on a tocopherol/tocotrienol scaffold",
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @classmethod
    def _is_tocopherol_like(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        return cls._carbon_count(participant) >= 28 and cls._oxygen_count(participant) == 2 and mol.GetRingInfo().NumRings() >= 2
