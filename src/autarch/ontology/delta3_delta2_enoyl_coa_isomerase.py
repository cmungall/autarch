"""delta(3)-delta(2)-enoyl-CoA isomerase activity.

Catalysis of migration of a double bond within an enoyl-CoA thioester from the
3-position to the 2-position, or the reverse direction.
"""

from collections import Counter

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.isomerase import Isomerase


class Delta3Delta2EnoylCoAIsomerase(Isomerase):
    """delta(3)-delta(2)-enoyl-CoA isomerase activity.

    Catalysis of migration of a double bond within an enoyl-CoA thioester from
    the 3-position to the 2-position, or the reverse direction.
    """

    GO_ID = "GO:0004165"
    EC_NUMBER_PREFIX = "5.3.3.8"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for reversible delta(3)/delta(2) double-bond migration in an enoyl-CoA."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not enoyl-CoA isomerase chemistry",
            )

        left_core = [participant for participant in reaction.left_participants if participant.get_mol() is not None]
        right_core = [participant for participant in reaction.right_participants if participant.get_mol() is not None]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected a single substrate/product pair for enoyl-CoA isomerization",
            )

        substrate = left_core[0]
        product = right_core[0]
        if not (self._is_coa_thioester(substrate) and self._is_coa_thioester(product)):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate and product are not both CoA-like thioesters",
            )
        if self._core_formula(substrate) != self._core_formula(product):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate and product do not preserve the enoyl-CoA scaffold",
            )

        if self._is_delta3_enoyl_coa(substrate) and self._is_delta2_enoyl_coa(product):
            return ClassificationResult(
                is_member=True,
                explanation="Enoyl-CoA isomerase: migration of a delta(3) double bond to the delta(2) position",
            )
        if self._is_delta2_enoyl_coa(substrate) and self._is_delta3_enoyl_coa(product):
            return ClassificationResult(
                is_member=True,
                explanation="Enoyl-CoA isomerase: reverse migration of a delta(2) double bond to the delta(3) position",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No delta(3)/delta(2) enoyl-CoA double-bond migration detected",
        )

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _is_coa_thioester(cls, participant: Participant) -> bool:
        return bool(participant.is_thioester() and cls._phosphorus_count(participant) >= 2)

    @staticmethod
    def _core_formula(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() not in {0, 1}
        )

    @classmethod
    def _acyl_chain_positions(cls, participant: Participant):
        mol = participant.get_mol()
        if mol is None or not cls._is_coa_thioester(participant):
            return None
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            has_thioester_sulfur = any(
                neighbor.GetAtomicNum() == 16
                and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 1
                for neighbor in atom.GetNeighbors()
            )
            has_carbonyl_oxygen = any(
                neighbor.GetAtomicNum() == 8
                and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2
                for neighbor in atom.GetNeighbors()
            )
            if not (has_thioester_sulfur and has_carbonyl_oxygen):
                continue
            alpha_candidates = [
                neighbor
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() == 6
            ]
            if len(alpha_candidates) != 1:
                continue
            alpha_carbon = alpha_candidates[0]
            beta_candidates = [
                neighbor
                for neighbor in alpha_carbon.GetNeighbors()
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != atom.GetIdx()
            ]
            if len(beta_candidates) != 1:
                continue
            beta_carbon = beta_candidates[0]
            gamma_candidates = [
                neighbor
                for neighbor in beta_carbon.GetNeighbors()
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_carbon.GetIdx()
            ]
            gamma_carbon = gamma_candidates[0] if len(gamma_candidates) == 1 else None
            return atom, alpha_carbon, beta_carbon, gamma_carbon
        return None

    @classmethod
    def _is_delta2_enoyl_coa(cls, participant: Participant) -> bool:
        positions = cls._acyl_chain_positions(participant)
        if positions is None:
            return False
        _, alpha_carbon, beta_carbon, _ = positions
        mol = participant.get_mol()
        if mol is None:
            return False
        alpha_beta_bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), beta_carbon.GetIdx())
        return alpha_beta_bond is not None and alpha_beta_bond.GetBondTypeAsDouble() == 2

    @classmethod
    def _is_delta3_enoyl_coa(cls, participant: Participant) -> bool:
        positions = cls._acyl_chain_positions(participant)
        if positions is None:
            return False
        _, alpha_carbon, beta_carbon, gamma_carbon = positions
        if gamma_carbon is None:
            return False
        if not cls._has_chain_extension(gamma_carbon, excluded_atom_ids={beta_carbon.GetIdx()}):
            return False
        mol = participant.get_mol()
        if mol is None:
            return False
        alpha_beta_bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), beta_carbon.GetIdx())
        beta_gamma_bond = mol.GetBondBetweenAtoms(beta_carbon.GetIdx(), gamma_carbon.GetIdx())
        return (
            alpha_beta_bond is not None
            and alpha_beta_bond.GetBondTypeAsDouble() == 1
            and beta_gamma_bond is not None
            and beta_gamma_bond.GetBondTypeAsDouble() == 2
        )

    @staticmethod
    def _has_chain_extension(atom, excluded_atom_ids: set[int]) -> bool:
        return any(
            neighbor.GetIdx() not in excluded_atom_ids and neighbor.GetAtomicNum() == 6
            for neighbor in atom.GetNeighbors()
        )
