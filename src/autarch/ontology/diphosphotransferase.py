"""diphosphotransferase.

Catalysis of the transfer of a diphosphate group from one compound (donor) to
another (acceptor).
"""

from __future__ import annotations

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_UDP,
    CHEBI_UTP,
)
from autarch.ontology.transferase import Transferase

CHEBI_GMP = "CHEBI:58115"
CHEBI_UMP = "CHEBI:57865"


class Diphosphotransferase(Transferase):
    """diphosphotransferase.

    Catalysis of the transfer of a diphosphate group from one compound (donor)
    to another (acceptor).
    """

    GO_ID = "GO:0016778"
    EC_NUMBER_PREFIX = "2.7.6.-"

    NTP_TO_NMP = {
        CHEBI_ATP: CHEBI_AMP,
        CHEBI_GTP: CHEBI_GMP,
        CHEBI_CTP: CHEBI_CMP,
        CHEBI_UTP: CHEBI_UMP,
    }
    NDP_PRODUCTS = {CHEBI_ADP, CHEBI_GDP, CHEBI_CDP, CHEBI_UDP}
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    P_ANHYDRIDE_PATTERN = Chem.MolFromSmarts("[P;!R]-O-[P;!R]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for NTP-driven diphosphate transfer to an acceptor."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not diphosphate transfer",
            )

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
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
        if any(participant.chebi_id == CHEBI_DIPHOSPHATE for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Free diphosphate release indicates ligase or donor-synthesis chemistry, not diphosphate transfer to an acceptor",
            )

        donor_pair = self._get_ntp_to_nmp_signature(left_participants, right_participants)
        if donor_pair is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing nucleoside triphosphate to nucleoside monophosphate donor signature",
            )
        donor_ntp, donor_nmp = donor_pair

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
            and participant.chebi_id != donor_ntp
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
            and participant.chebi_id != donor_nmp
            and participant.chebi_id not in self.NDP_PRODUCTS
        ]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No acceptor/product pair remained after removing NTP/NMP donors and solvent spectators",
            )

        for substrate in left_reactive:
            for product in right_reactive:
                if self._is_diphosphorylated_pair(substrate, product):
                    return ClassificationResult(
                        is_member=True,
                        explanation="Diphosphotransferase: an NTP donor transfers a diphosphate group to an acceptor with NMP release",
                    )

        return ClassificationResult(
            is_member=False,
            explanation="No acceptor/product pair showed diphosphate transfer chemistry",
        )

    def _get_ntp_to_nmp_signature(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> tuple[str, str] | None:
        for ntp, nmp in self.NTP_TO_NMP.items():
            has_ntp = any(participant.chebi_id == ntp for participant in left_participants)
            has_nmp = any(participant.chebi_id == nmp for participant in right_participants)
            if has_ntp and has_nmp:
                return ntp, nmp
        return None

    def _is_diphosphorylated_pair(
        self,
        substrate: Participant,
        product: Participant,
    ) -> bool:
        if substrate.get_mol() is None or product.get_mol() is None:
            return False
        if self._non_phosphorus_scaffold_counts(substrate) != self._non_phosphorus_scaffold_counts(product):
            return False
        if self._phosphorus_count(product) <= self._phosphorus_count(substrate):
            return False
        return self._has_p_anhydride(product) and not self._has_p_anhydride(substrate)

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @staticmethod
    def _non_phosphorus_scaffold_counts(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() not in {1, 8, 15}
        )

    @classmethod
    def _has_p_anhydride(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        return mol is not None and cls.P_ANHYDRIDE_PATTERN is not None and mol.HasSubstructMatch(cls.P_ANHYDRIDE_PATTERN)
