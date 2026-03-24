"""aminoacyltransferase.

Catalysis of the transfer of an amino-acyl group from one compound (donor) to
 another (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    )
from autarch.ontology.transferase import Transferase


class Aminoacyltransferase(Transferase):
    """aminoacyltransferase.

    Catalysis of the transfer of an amino-acyl group from one compound (donor)
    to another (acceptor).
    """

    GO_ID = "GO:0016755"
    EC_NUMBER_PREFIX = "2.3.2.-"
    ACYL_ESTER_PATTERN = Chem.MolFromSmarts("[O;X2][C;X3](=[OX1])[C,N]")
    AMIDE_PATTERN = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    GLUTAMINE_DONOR_PATTERN = Chem.MolFromSmarts("[NX3H2,X3H1][CX3](=[OX1])[CH2][CH2][CH]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for aminoacyl transfer reactions."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not aminoacyl transfer",
            )

        all_chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if CHEBI_ATP in all_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="ATP-dependent bond formation indicates ligase chemistry rather than aminoacyl transfer",
            )
        if CHEBI_O2 in all_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen chemistry indicates an oxidoreductase rather than aminoacyl transferase",
            )
        if {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS, CHEBI_NADH, CHEBI_NADPH} & all_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="Nicotinamide cofactors indicate redox chemistry rather than aminoacyl transfer",
            )

        trna_branch = self._has_aminoacyl_trna_transfer(reaction)
        if trna_branch:
            return ClassificationResult(
                is_member=True,
                explanation="Aminoacyltransferase: aminoacyl group transferred from aminoacyl-tRNA to an acceptor",
            )

        amide_branch = self._has_amide_transfer_branch(reaction)
        if amide_branch:
            return ClassificationResult(
                is_member=True,
                explanation="Aminoacyltransferase: aminoacyl group transferred from an amide donor with ammonium release",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No aminoacyl-transfer signature detected",
        )

    def _has_aminoacyl_trna_transfer(self, reaction: Reaction) -> bool:
        """Detect aminoacyl-tRNA donor chemistry."""
        left_trna_donors = [
            participant
            for participant in reaction.left_participants
            if self._is_trna_like(participant) and self._has_pattern(participant, self.ACYL_ESTER_PATTERN)
        ]
        if not left_trna_donors:
            return False

        right_trna_products = [
            participant
            for participant in reaction.right_participants
            if self._is_trna_like(participant)
        ]
        if not right_trna_products:
            return False

        left_acceptors = [
            participant
            for participant in reaction.left_participants
            if not self._is_trna_like(participant) and self._is_substantive_carbon_acceptor(participant)
        ]
        right_acceptors = [
            participant
            for participant in reaction.right_participants
            if not self._is_trna_like(participant)
            and participant.chebi_id != CHEBI_H_PLUS
            and self._is_substantive_carbon_acceptor(participant)
        ]
        return bool(left_acceptors) and bool(right_acceptors)

    def _has_amide_transfer_branch(self, reaction: Reaction) -> bool:
        """Detect glutamyl/glutaminyl transfer with ammonium release."""
        if any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return False
        if not any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in reaction.right_participants):
            return False

        left_reactive = [participant for participant in reaction.left_participants if participant.chebi_id != CHEBI_H_PLUS]
        right_reactive = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in {CHEBI_H_PLUS, CHEBI_NH3, CHEBI_NH4}
        ]
        if not right_reactive:
            return False

        if left_reactive and right_reactive and all(
            self._is_protein_like(participant) for participant in left_reactive + right_reactive
        ):
            return True

        amide_donors = [
            participant
            for participant in left_reactive
            if self._has_pattern(participant, self.GLUTAMINE_DONOR_PATTERN)
        ]
        if not amide_donors:
            return False

        if len(left_reactive) < 2:
            return False

        return any(self._has_pattern(participant, self.AMIDE_PATTERN) for participant in right_reactive)

    @staticmethod
    def _is_trna_like(participant: Participant) -> bool:
        return participant.polymer_type == PolymerType.TRNA or participant.chebi_id == "CHEBI:78442"

    @staticmethod
    def _is_protein_like(participant: Participant) -> bool:
        return participant.polymer_type in {PolymerType.PROTEIN, PolymerType.PEPTIDE}

    @staticmethod
    def _has_pattern(participant: Participant, pattern: Chem.Mol | None) -> bool:
        if pattern is None:
            return False
        mol = participant.get_mol()
        return mol is not None and mol.HasSubstructMatch(pattern)

    @staticmethod
    def _is_substantive_carbon_acceptor(participant: Participant) -> bool:
        if participant.polymer_type in {PolymerType.PROTEIN, PolymerType.PEPTIDE}:
            return True
        mol = participant.get_mol()
        if mol is None:
            return False
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) > 0
