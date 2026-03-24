"""L-histidine N(alpha)-methyltransferase activity.

Catalysis of transfer of a methyl group from S-adenosyl-L-methionine to the
alpha-amino nitrogen of a histidine-derived amino acid substrate.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase


class LHistidineNAlphaMethyltransferase(Methyltransferase):
    """L-histidine N(alpha)-methyltransferase activity.

    Catalysis of transfer of a methyl group from S-adenosyl-L-methionine to the
    alpha-amino nitrogen of a histidine-derived amino acid substrate.
    """

    GO_ID = "GO:0052706"
    EC_NUMBER_PREFIX = "2.1.1.44"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for SAM-dependent methylation of a histidine-like amino acid."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}
        ]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one histidine-derived substrate/product pair",
            )

        substrate = left_core[0]
        product = right_core[0]
        if not self._is_histidine_like_amino_acid(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="No histidine-like amino acid substrate detected",
            )
        if not self._is_histidine_like_amino_acid(product):
            return ClassificationResult(
                is_member=False,
                explanation="No histidine-like amino acid product detected",
            )
        if self._carbon_count(product) != self._carbon_count(substrate) + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Product is not a single-methyl homolog of the histidine-like substrate",
            )
        if self._nitrogen_count(product) != self._nitrogen_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Histidine nitrogen count is not preserved across methyl transfer",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Histidine N(alpha)-methyltransferase: SAM-dependent methylation of a histidine-derived amino acid",
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @classmethod
    def _is_histidine_like_amino_acid(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        if not participant.has_moiety(Moiety.CARBOXYL):
            return False
        if cls._oxygen_count(participant) != 2:
            return False
        if cls._nitrogen_count(participant) < 3:
            return False
        if cls._carbon_count(participant) > 12:
            return False
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) != 5:
                continue
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if sum(atom.GetAtomicNum() == 7 for atom in atoms) >= 2:
                return True
        return False
