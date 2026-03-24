"""Catalysis of the hydrolysis of an ether bond, -O-.

This EC 3.3.2 branch covers hydrolytic cleavage of non-phosphoric ether bonds,
including epoxide and oxetane ring opening as well as ether cleavage that yields
more oxidized oxygenated products.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_CTP,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.hydrolase import Hydrolase


class EtherHydrolase(Hydrolase):
    """ether hydrolase

    Catalysis of the hydrolysis of an ether bond, -O-.

    The class is defined by water-dependent cleavage of an ether-like linkage,
    especially epoxide or oxetane ring opening and related ether hydrolysis that
    increases alcohol or carbonyl functionality in the products.
    """

    GO_ID = "GO:0016803"
    EC_NUMBER_PREFIX = "3.3.2.-"

    EXCLUDED_CHEBIS = {
        CHEBI_ATP,
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_GTP,
        CHEBI_GDP,
        CHEBI_CTP,
        CHEBI_CDP,
        CHEBI_CMP,
        CHEBI_UTP,
        CHEBI_UDP,
        CHEBI_UMP,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H_PLUS}
    ETHER_PATTERN = Chem.MolFromSmarts("[C;!$(C=O)][O;X2][C;!$(C=O)]")
    EPOXIDE_PATTERN = Chem.MolFromSmarts("[O;r3]1[C;r3][C;r3]1")
    OXETANE_PATTERN = Chem.MolFromSmarts("[O;r4]1[C;r4][C;r4][C;r4]1")
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1][C;!$(C=O)]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")
    ALDEHYDE_PATTERN = Chem.MolFromSmarts("[CX3H1,H2](=O)")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolytic ether cleavage."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not ether hydrolysis",
            )

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External redox or ATP coupling indicates another reaction class",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Water is required as the hydrolytic cosubstrate",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No chemically resolved ether substrate/product pair remains after removing solvent spectators",
            )

        if self._has_small_ring_ether_opening(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Ether hydrolase: hydrolytic opening of an epoxide or oxetane ring",
            )

        if self._has_general_ether_cleavage(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Ether hydrolase: hydrolysis decreases ether bonding and increases product oxygenation",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No hydrolytic ether-cleavage signature detected",
        )

    def _has_small_ring_ether_opening(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        left_ring_ethers = sum(
            self._substructure_count(participant, self.EPOXIDE_PATTERN)
            + self._substructure_count(participant, self.OXETANE_PATTERN)
            for participant in left_reactive
        )
        right_ring_ethers = sum(
            self._substructure_count(participant, self.EPOXIDE_PATTERN)
            + self._substructure_count(participant, self.OXETANE_PATTERN)
            for participant in right_reactive
        )
        left_hydroxyls = sum(self._substructure_count(participant, self.HYDROXYL_PATTERN) for participant in left_reactive)
        right_hydroxyls = sum(self._substructure_count(participant, self.HYDROXYL_PATTERN) for participant in right_reactive)
        return left_ring_ethers > right_ring_ethers and right_hydroxyls > left_hydroxyls

    def _has_general_ether_cleavage(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        left_ethers = sum(self._substructure_count(participant, self.ETHER_PATTERN) for participant in left_reactive)
        right_ethers = sum(self._substructure_count(participant, self.ETHER_PATTERN) for participant in right_reactive)
        left_carbonyls = sum(self._substructure_count(participant, self.CARBONYL_PATTERN) for participant in left_reactive)
        right_carbonyls = sum(self._substructure_count(participant, self.CARBONYL_PATTERN) for participant in right_reactive)
        if not (left_ethers > right_ethers and right_carbonyls > left_carbonyls):
            return False

        right_aldehydes = sum(self._substructure_count(participant, self.ALDEHYDE_PATTERN) for participant in right_reactive)
        left_aldehydes = sum(self._substructure_count(participant, self.ALDEHYDE_PATTERN) for participant in left_reactive)
        return right_aldehydes > left_aldehydes or any(
            0 < self._carbon_count(participant) <= 3 for participant in right_reactive
        )

    @staticmethod
    def _substructure_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) * max(participant.count, 1)
