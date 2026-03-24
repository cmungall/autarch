"""Catalysis of the transfer of an acyl group, other than amino-acyl, from one compound (donor) to another (acceptor)."""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_COA,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.transferase import Transferase


class AcyltransferaseTransferringGroupsOtherThanAminoAcylGroups(Transferase):
    """acyltransferase transferring groups other than amino-acyl groups

    Catalysis of the transfer of an acyl group, other than amino-acyl, from one compound (donor) to another (acceptor).
    """

    GO_ID = "GO:0016747"
    EC_NUMBER_PREFIX = "2.3.1.-"
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H_PLUS}
    THIOESTER_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][PX4](=[OX1])([OX2,OX1-])[OX2,OX1-]")
    ESTER_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")
    SUGAR_RING_PATTERN = Chem.MolFromSmarts("[O;R]1[CH;R][CH;R][CH;R][CH;R][CH;R]1")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for non-aminoacyl acyl transfer."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not acyl transfer",
            )

        chebi_ids = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if self._has_triphosphate_coupling(reaction.left_participants, reaction.right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="ATP/GTP-coupled bond formation indicates ligase chemistry",
            )
        if CHEBI_O2 in chebi_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen chemistry indicates a different enzyme class",
            )
        if self._has_redox_pair(chebi_ids):
            transferable_acyl_count = sum(
                1
                for participant in reaction.left_participants + reaction.right_participants
                if self._is_transferable_acyl_donor(participant)
            )
            if transferable_acyl_count < 2:
                return ClassificationResult(
                    is_member=False,
                    explanation="NAD(P)-linked acyl activation without an acyl-transfer scaffold indicates a different class",
                )

        forward = self._check_direction(reaction.left_participants, reaction.right_participants)
        if forward.is_member:
            return forward

        reverse = self._check_direction(reaction.right_participants, reaction.left_participants)
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        left = [participant for participant in left_participants if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS]
        right = [participant for participant in right_participants if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS]
        if any(participant.chebi_id == CHEBI_H2O for participant in left_participants) and len(left) == 1:
            return ClassificationResult(
                is_member=False,
                explanation="Simple hydrolysis is not acyl transfer",
            )

        if CHEBI_SAM in {participant.chebi_id for participant in left_participants} and CHEBI_SAH in {
            participant.chebi_id for participant in right_participants
        }:
            return ClassificationResult(
                is_member=False,
                explanation="SAM-dependent alkyl transfer is not acyltransferase chemistry",
            )

        has_acyl_donor = any(self._is_transferable_acyl_donor(participant) for participant in left)
        if not has_acyl_donor:
            return ClassificationResult(
                is_member=False,
                explanation="No transferable acyl donor detected",
            )

        if self._is_coa_transferase_like(left, right):
            return ClassificationResult(
                is_member=False,
                explanation="Acyl carrier remains thioester-bound on the product side, indicating CoA-transferase-like chemistry",
            )

        if self._is_water_assisted_thioester_condensation(left_participants, left, right):
            return ClassificationResult(
                is_member=False,
                explanation="Water-assisted thioester condensation indicates a different acyl-transfer branch",
            )

        if self._has_carrier_release(right) and self._has_new_acylated_product(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Acyltransferase: donor acyl group is transferred to a new acceptor with carrier release",
            )

        if self._is_decarboxylative_condensation(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Acyltransferase: acyl donor undergoes carrier release during decarboxylative acyl transfer",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No non-aminoacyl acyl-transfer signature detected",
        )

    def _is_coa_transferase_like(self, left: list[Participant], right: list[Participant]) -> bool:
        """Reject reactions that move an intact thioester carrier between acids."""
        return (
            any(participant.is_thioester() for participant in left)
            and any(participant.is_thioester() for participant in right)
            and not self._has_carrier_release(right)
        )

    def _is_water_assisted_thioester_condensation(
        self,
        left_participants: list[Participant],
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Reject water-assisted C-C condensations that form a new thioester product."""
        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return False
        if any(participant.chebi_id == CHEBI_CO2 for participant in right):
            return False
        new_products = [
            participant
            for participant in right
            if self._is_acylated_product(participant) and not any(participant.is_same_molecule(other) for other in left)
        ]
        return bool(new_products) and all(participant.is_thioester() for participant in new_products)

    def _has_new_acylated_product(self, left: list[Participant], right: list[Participant]) -> bool:
        """Detect amide, ester, or thioester formation on a new product scaffold."""
        return any(
            self._is_acylated_product(participant) and not any(participant.is_same_molecule(other) for other in left)
            for participant in right
        )

    def _is_decarboxylative_condensation(self, left: list[Participant], right: list[Participant]) -> bool:
        """Detect carrier-releasing condensations within EC 2.3.1."""
        return (
            len(left) >= 2
            and any(participant.is_thioester() or self._matches(participant, self.ACYL_PHOSPHATE_PATTERN) for participant in left)
            and self._has_carrier_release(right)
            and any(participant.chebi_id == CHEBI_CO2 for participant in right)
            and any(
                self._matches(participant, self.CARBONYL_PATTERN)
                and not any(participant.is_same_molecule(other) for other in left)
                for participant in right
            )
        )

    def _has_carrier_release(self, participants: list[Participant]) -> bool:
        """Detect release of CoA, phosphate, sulfur carriers, or sugar carriers."""
        return any(
            participant.chebi_id in {CHEBI_COA, CHEBI_PHOSPHATE}
            or self._is_released_sulfur_carrier(participant)
            or self._is_sugar_like_release(participant)
            for participant in participants
        )

    def _is_transferable_acyl_donor(self, participant: Participant) -> bool:
        """Detect thioester, acyl-phosphate, or activated ester donors."""
        return (
            participant.is_thioester()
            or self._matches(participant, self.ACYL_PHOSPHATE_PATTERN)
            or self._is_activated_ester_donor(participant)
        )

    def _is_activated_ester_donor(self, participant: Participant) -> bool:
        """Detect ester donors such as acyl-sugars."""
        return self._matches(participant, self.ESTER_PATTERN) and self._sugar_ring_count(participant) >= 1

    def _is_acylated_product(self, participant: Participant) -> bool:
        """Detect products carrying the transferred acyl group."""
        return (
            participant.is_thioester()
            or self._matches(participant, self.AMIDE_PATTERN)
            or self._matches(participant, self.ESTER_PATTERN)
        )

    @staticmethod
    def _matches(participant: Participant, pattern: Chem.Mol | None) -> bool:
        """Check whether a participant matches a SMARTS pattern."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)

    @staticmethod
    def _has_triphosphate_coupling(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect ATP- or GTP-dependent ligase chemistry."""
        left_ids = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        return (
            CHEBI_ATP in left_ids and bool({CHEBI_ADP, CHEBI_AMP} & right_ids)
        ) or (
            CHEBI_GTP in left_ids and CHEBI_GDP in right_ids
        )

    @staticmethod
    def _is_released_sulfur_carrier(participant: Participant) -> bool:
        """Detect sulfur-containing carrier release without a remaining thioester bond."""
        mol = participant.get_mol()
        if mol is None or participant.is_thioester():
            return False
        return any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())

    def _is_sugar_like_release(self, participant: Participant) -> bool:
        """Detect sugar leaving groups such as glucose released from acyl-sugars."""
        return (
            self._sugar_ring_count(participant) >= 1
            and self._phosphorus_count(participant) == 0
            and self._nitrogen_count(participant) == 0
        )

    def _sugar_ring_count(self, participant: Participant) -> int:
        """Count hexose-like sugar rings."""
        mol = participant.get_mol()
        if mol is None or self.SUGAR_RING_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(self.SUGAR_RING_PATTERN))

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        """Count phosphorus atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @staticmethod
    def _has_redox_pair(chebi_ids: set[str]) -> bool:
        """Detect NAD(P)-linked redox chemistry."""
        return (
            CHEBI_NAD_PLUS in chebi_ids and CHEBI_NADH in chebi_ids
        ) or (
            CHEBI_NADP_PLUS in chebi_ids and CHEBI_NADPH in chebi_ids
        )

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        """Count nitrogen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
