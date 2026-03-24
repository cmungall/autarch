"""carboxyl or carbamoyltransferase.

Catalysis of the transfer of a carboxyl- or carbamoyl group from one compound
(donor) to another (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_PHOSPHATE
from autarch.moiety import is_ketone
from autarch.ontology.transferase import Transferase


class CarboxylOrCarbamoyltransferase(Transferase):
    """carboxyl or carbamoyltransferase.

    Catalysis of the transfer of a carboxyl- or carbamoyl group from one
    compound (donor) to another (acceptor).
    """

    GO_ID = "GO:0016743"
    EC_NUMBER_PREFIX = "2.1.3.-"
    CARBAMOYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[OX2][PX4]")
    AMINE_PATTERN = Chem.MolFromSmarts("[N;H3+,H2,H1;!$(N=*);!$(N#*)]")
    CARBAMOYL_PRODUCT_PATTERN = Chem.MolFromSmarts("[NX3,O][CX3](=[OX1])[NX3]")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,-]")
    PHOSPHATE_PRODUCT_CHEBIS = {CHEBI_PHOSPHATE, "CHEBI:16838"}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for carboxyl- or carbamoyl-group transfer."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not carboxyl or carbamoyl transfer",
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
        if self._is_carbamoyl_transfer(left_participants, right_participants):
            return ClassificationResult(
                is_member=True,
                explanation="Carboxyl or carbamoyltransferase: carbamoyl phosphate donates a carbamoyl group to an acceptor substrate",
            )

        if self._is_carboxyl_transfer(left_participants, right_participants):
            return ClassificationResult(
                is_member=True,
                explanation="Carboxyl or carbamoyltransferase: a carboxyl group shifts from an activated donor to an acceptor substrate",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No carboxyl- or carbamoyl-transfer signature detected",
        )

    def _is_carbamoyl_transfer(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        left_donors = [
            participant
            for participant in left_participants
            if self._matches(participant, self.CARBAMOYL_PHOSPHATE_PATTERN)
        ]
        if not left_donors:
            return False

        if not any(participant.chebi_id in self.PHOSPHATE_PRODUCT_CHEBIS for participant in right_participants):
            return False

        left_acceptors = [
            participant
            for participant in left_participants
            if participant not in left_donors and participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_products = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.PHOSPHATE_PRODUCT_CHEBIS | {CHEBI_H_PLUS}
        ]
        if not left_acceptors or not right_products:
            return False

        if not any(self._matches(participant, self.AMINE_PATTERN) for participant in left_acceptors):
            return False

        return self._pattern_total(right_products, self.CARBAMOYL_PRODUCT_PATTERN) > self._pattern_total(
            left_acceptors, self.CARBAMOYL_PRODUCT_PATTERN
        )

    def _is_carboxyl_transfer(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        left = [participant for participant in left_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        right = [participant for participant in right_participants if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}]
        left_thioesters = [participant for participant in left if participant.is_thioester()]
        right_thioesters = [participant for participant in right if participant.is_thioester()]
        left_acceptors = [participant for participant in left if self._is_keto_acid(participant)]
        right_acceptors = [participant for participant in right if self._is_keto_acid(participant)]
        if not left_thioesters or not right_thioesters or not left_acceptors or not right_acceptors:
            return False

        if self._carboxyl_total(left) != self._carboxyl_total(right):
            return False

        donor_loses_carboxyl = max(self._carboxyl_count(participant) for participant in left_thioesters) > max(
            self._carboxyl_count(participant) for participant in right_thioesters
        )
        acceptor_gains_carboxyl = max(self._carboxyl_count(participant) for participant in right_acceptors) > max(
            self._carboxyl_count(participant) for participant in left_acceptors
        )
        return donor_loses_carboxyl and acceptor_gains_carboxyl

    @staticmethod
    def _matches(participant: Participant, pattern: Chem.Mol | None) -> bool:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)

    @classmethod
    def _pattern_total(cls, participants: list[Participant], pattern: Chem.Mol | None) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None or pattern is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total

    def _carboxyl_total(self, participants: list[Participant]) -> int:
        return sum(self._carboxyl_count(participant) * max(participant.count, 1) for participant in participants)

    def _carboxyl_count(self, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return len(mol.GetSubstructMatches(self.CARBOXYL_PATTERN))

    @staticmethod
    def _is_keto_acid(participant: Participant) -> bool:
        return bool(
            participant.smiles
            and is_ketone(participant.smiles, participant.chebi_id)
            and participant.get_mol() is not None
        )
