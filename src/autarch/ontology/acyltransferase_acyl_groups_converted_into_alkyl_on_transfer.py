"""acyltransferase acyl groups converted into alkyl on transfer.

Catalysis of the transfer of an acyl group from one compound (donor) to another
(acceptor), with the acyl group being converted into alkyl on transfer.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_H2O2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.moiety import is_aldehyde, is_carboxylic_acid, is_ketone
from autarch.ontology.transferase import Transferase


class AcyltransferaseAcylGroupsConvertedIntoAlkylOnTransfer(Transferase):
    """acyltransferase acyl groups converted into alkyl on transfer.

    Catalysis of the transfer of an acyl group from one compound (donor) to
    another (acceptor), with the acyl group being converted into alkyl on
    transfer.
    """

    GO_ID = "GO:0046912"
    EC_NUMBER_PREFIX = "2.3.3.-"
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts(
        "[CX3](=[OX1])[OX2][PX4](=[OX1])([OX2,OX1-])[OX2,OX1-]"
    )
    AMIDE_PATTERN = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for C-C acyl transfer with carrier release."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not acyl transfer",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & {
            CHEBI_ATP,
            CHEBI_ADP,
            CHEBI_AMP,
            CHEBI_GTP,
            CHEBI_GDP,
            CHEBI_NAD_PLUS,
            CHEBI_NADH,
            CHEBI_NADP_PLUS,
            CHEBI_NADPH,
            CHEBI_FAD,
            CHEBI_FADH2,
            CHEBI_O2,
            CHEBI_H2O2,
        }:
            return ClassificationResult(
                is_member=False,
                explanation="Energy-coupled, redox, or oxygen chemistry indicates another transferase branch",
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
        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Condensing acyltransferase chemistry requires water in this operational approximation",
            )

        left = [participant for participant in left_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        right = [participant for participant in right_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        donors = [participant for participant in left if self._is_activated_acyl_donor(participant)]
        if not donors:
            return ClassificationResult(
                is_member=False,
                explanation="No activated acyl donor detected",
            )

        if not self._has_carrier_release(right):
            return ClassificationResult(
                is_member=False,
                explanation="No carrier-release product detected",
            )

        acceptors = [
            participant
            for participant in left
            if participant not in donors and self._is_condensation_acceptor(participant)
        ]
        if not acceptors:
            acceptors = [participant for participant in left if self._is_condensation_acceptor(participant)]
        if not acceptors:
            return ClassificationResult(
                is_member=False,
                explanation="No carbonyl-bearing acceptor substrate detected",
            )

        new_products = [
            participant
            for participant in right
            if not self._has_carrier_release([participant])
            and not any(participant.is_same_molecule(other) for other in left)
        ]
        if not new_products:
            return ClassificationResult(
                is_member=False,
                explanation="No new condensation product detected",
            )

        max_acceptor_c = max(self._carbon_count(participant) for participant in acceptors)
        if not any(
            self._carbon_count(participant) > max_acceptor_c
            and (
                participant.is_thioester()
                or self._is_acid_like(participant)
                or self._matches(participant, self.AMIDE_PATTERN)
            )
            for participant in new_products
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Products do not show the larger acyl-condensation scaffold expected for EC 2.3.3 chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Acyltransferase: an activated acyl donor undergoes condensation with a carbonyl acceptor, converting the transferred acyl unit into alkyl on transfer",
        )

    def _is_activated_acyl_donor(self, participant: Participant) -> bool:
        return participant.is_thioester() or self._matches(participant, self.ACYL_PHOSPHATE_PATTERN)

    def _is_condensation_acceptor(self, participant: Participant) -> bool:
        if participant.is_thioester():
            return True
        if not participant.smiles:
            return False
        return is_aldehyde(participant.smiles, participant.chebi_id) or is_ketone(
            participant.smiles, participant.chebi_id
        )

    def _is_acid_like(self, participant: Participant) -> bool:
        return bool(participant.smiles and is_carboxylic_acid(participant.smiles, participant.chebi_id))

    @staticmethod
    def _matches(participant: Participant, pattern: Chem.Mol | None) -> bool:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _has_carrier_release(participants: list[Participant]) -> bool:
        for participant in participants:
            if participant.chebi_id == CHEBI_PHOSPHATE:
                return True
            mol = participant.get_mol()
            if mol is None or participant.is_thioester():
                continue
            if any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
                return True
        return False
