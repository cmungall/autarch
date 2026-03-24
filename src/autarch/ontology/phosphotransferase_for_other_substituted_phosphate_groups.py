"""phosphotransferase activity, for other substituted phosphate groups.

Catalysis of the transfer of a substituted phosphate-containing group from one compound to an acceptor.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_GDP,
    CHEBI_GMP,
    CHEBI_GTP,
    CHEBI_H_PLUS,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.transferase import Transferase

CHEBI_PAP = "CHEBI:58343"


class PhosphotransferaseForOtherSubstitutedPhosphateGroups(Transferase):
    """phosphotransferase activity, for other substituted phosphate groups.

    Catalysis of the transfer of a substituted phosphate-containing group from
    one compound to an acceptor.
    """

    GO_ID = "GO:0016780"
    EC_NUMBER_PREFIX = "2.7.8.-"
    LEAVING_GROUPS = {CHEBI_CMP, CHEBI_UMP, CHEBI_GMP, CHEBI_PAP}
    SIMPLE_PHOSPHATE_DONORS = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP}
    SIMPLE_SPENT_DONORS = {CHEBI_ADP, CHEBI_GDP, CHEBI_UDP, CHEBI_CDP}
    EXCLUDED_CHEBIS = SIMPLE_PHOSPHATE_DONORS | SIMPLE_SPENT_DONORS | {CHEBI_H_PLUS}
    FALLBACK_PHOSPHORUS_COUNTS = {
        CHEBI_CMP: 1,
        CHEBI_UMP: 1,
        CHEBI_GMP: 1,
        CHEBI_PAP: 2,
        CHEBI_ATP: 3,
        CHEBI_ADP: 2,
        CHEBI_GTP: 3,
        CHEBI_GDP: 2,
        CHEBI_UTP: 3,
        CHEBI_UDP: 2,
        CHEBI_CDP: 2,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for transfer of a substituted phosphate group."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not substituted-phosphate transfer",
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
        """Evaluate one transfer direction."""
        right_ids = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if right_ids & self.SIMPLE_SPENT_DONORS:
            return ClassificationResult(
                is_member=False,
                explanation="Simple nucleoside diphosphate byproducts indicate a different phosphotransferase branch",
            )

        leaving_group = next(
            (
                participant
                for participant in right_participants
                if participant.chebi_id in self.LEAVING_GROUPS
            ),
            None,
        )
        if leaving_group is None:
            return ClassificationResult(
                is_member=False,
                explanation="No substituted-phosphate leaving group detected",
            )

        donor_candidates = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.EXCLUDED_CHEBIS
            and self._phosphorus_count(participant) >= self._phosphorus_count(leaving_group)
        ]
        if not donor_candidates:
            return ClassificationResult(
                is_member=False,
                explanation="No substituted-phosphate donor candidate detected",
            )

        for donor in donor_candidates:
            left_reactive = [
                participant
                for participant in left_participants
                if participant is not donor and participant.chebi_id != CHEBI_H_PLUS
            ]
            right_reactive = [
                participant
                for participant in right_participants
                if participant is not leaving_group and participant.chebi_id != CHEBI_H_PLUS
            ]
            if not left_reactive or not right_reactive:
                continue
            left_total_p = sum(self._phosphorus_count(participant) for participant in left_reactive)
            right_total_p = sum(self._phosphorus_count(participant) for participant in right_reactive)
            donor_p = self._phosphorus_count(donor)
            leaving_group_p = self._phosphorus_count(leaving_group)
            if right_total_p != left_total_p + donor_p - leaving_group_p:
                continue
            if right_total_p <= left_total_p:
                continue
            if not any(self._phosphorus_count(participant) > 0 for participant in right_reactive):
                continue
            return ClassificationResult(
                is_member=True,
                explanation=(
                    "Phosphotransferase for other substituted phosphate groups: "
                    "a phosphorus-rich donor transfers a substituted phosphate "
                    "moiety while releasing a monophosphate nucleotide or PAP"
                ),
            )

        return ClassificationResult(
            is_member=False,
            explanation="No conserved substituted-phosphate transfer pattern detected",
        )

    def _phosphorus_count(self, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is not None:
            return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15) * max(
                participant.count,
                1,
            )
        if participant.chebi_id in self.FALLBACK_PHOSPHORUS_COUNTS:
            return self.FALLBACK_PHOSPHORUS_COUNTS[participant.chebi_id] * max(
                participant.count,
                1,
            )
        return 0
