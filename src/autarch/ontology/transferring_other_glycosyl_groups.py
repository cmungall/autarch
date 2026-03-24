"""transferring other glycosyl groups.

Operational EC-only approximation for EC 2.4.99 chemistry. This class is kept
precision-biased because the branch mixes Kdo transfer, heptose transfer, and
lipid-linked oligosaccharide transfer, and some RHEA records lack full
structures for the transferred glycosyl group.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_CMP,
    CHEBI_GDP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
)
from autarch.ontology.reaction import ReactionClass

CHEBI_DOLICHYL_DIPHOSPHATE = "CHEBI:57497"
CHEBI_DOLICHYL_DIPHOSPHOOLIGOSACCHARIDE = "CHEBI:57570"
CHEBI_DOLICHYL_GLYCAN = "CHEBI:132523"


class TransferringOtherGlycosylGroups(ReactionClass):
    """transferring other glycosyl groups.

    Operational EC-only approximation for EC 2.4.99 chemistry. This class is
    kept precision-biased because the branch mixes Kdo transfer, heptose
    transfer, and lipid-linked oligosaccharide transfer, and some RHEA records
    lack full structures for the transferred glycosyl group.
    """

    GO_ID = None
    EC_NUMBER_PREFIX = "2.4.99.-"
    SPECTATOR_IDS = {CHEBI_H2O, CHEBI_H_PLUS}
    CARRIER_CARBONS = {
        CHEBI_CMP: 9,
        CHEBI_ADP: 10,
        CHEBI_GDP: 10,
    }
    DOLICHYL_OLIGO_DONOR_IDS = {
        CHEBI_DOLICHYL_DIPHOSPHOOLIGOSACCHARIDE,
        CHEBI_DOLICHYL_GLYCAN,
    }
    DOLICHYL_CARRIER_IDS = {CHEBI_DOLICHYL_DIPHOSPHATE}
    RING_OXYGEN_PATTERN = Chem.MolFromSmarts("[O;R]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for EC 2.4.99 glycosyl transfer chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not glycosyl transfer",
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one transfer orientation."""
        if self._is_lipid_linked_oligosaccharide_transfer(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Transferring other glycosyl groups: lipid-linked oligosaccharide transfer to a protein acceptor",
            )

        if self._is_other_nucleotide_glycosyl_transfer(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Transferring other glycosyl groups: a non-pentose, non-hexose glycosyl unit is transferred from a carrier-bound donor",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No EC 2.4.99 glycosyl-transfer signature detected",
        )

    def _is_other_nucleotide_glycosyl_transfer(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Detect CMP-Kdo, ADP-heptose, and similar donor chemistries."""
        right_carriers = [
            participant for participant in right if self._carrier_carbons(participant) is not None
        ]
        if not right_carriers:
            return False

        left_reactive = [
            participant
            for participant in left
            if participant.chebi_id not in self.SPECTATOR_IDS
        ]
        right_reactive = [
            participant
            for participant in right
            if participant.chebi_id not in self.SPECTATOR_IDS
        ]

        for carrier in right_carriers:
            carrier_carbons = self._carrier_carbons(carrier)
            if carrier_carbons is None:
                continue
            carrier_phosphorus = self._phosphorus_count(carrier)
            for donor in left_reactive:
                donor_mol = donor.get_mol()
                if donor_mol is None:
                    continue
                donor_carbons = self._carbon_count(donor)
                donor_phosphorus = self._phosphorus_count(donor)
                donor_transfer_carbons = donor_carbons - carrier_carbons
                if donor_transfer_carbons not in {7, 8}:
                    continue
                if carrier.chebi_id == CHEBI_CMP and donor_transfer_carbons != 8:
                    continue
                if carrier.chebi_id in {CHEBI_ADP, CHEBI_GDP} and donor_transfer_carbons != 7:
                    continue
                if donor_phosphorus < carrier_phosphorus:
                    continue
                carrier_nitrogen = self._nitrogen_count(carrier)
                if self._nitrogen_count(donor) != carrier_nitrogen:
                    continue
                for acceptor in left_reactive:
                    if acceptor is donor or acceptor.get_mol() is None:
                        continue
                    acceptor_carbons = self._carbon_count(acceptor)
                    acceptor_phosphorus = self._phosphorus_count(acceptor)
                    acceptor_ring_oxygens = self._ring_oxygen_count(acceptor)
                    for product in right_reactive:
                        if product is carrier or product.get_mol() is None:
                            continue
                        if self._carbon_count(product) != acceptor_carbons + donor_transfer_carbons:
                            continue
                        if self._phosphorus_count(product) != (
                            acceptor_phosphorus + donor_phosphorus - carrier_phosphorus
                        ):
                            continue
                        if self._ring_oxygen_count(product) <= acceptor_ring_oxygens:
                            continue
                        return True
        return False

    def _is_lipid_linked_oligosaccharide_transfer(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Detect dolichyl-linked oligosaccharide transfer when structures are absent."""
        has_donor = any(participant.chebi_id in self.DOLICHYL_OLIGO_DONOR_IDS for participant in left)
        has_carrier = any(participant.chebi_id in self.DOLICHYL_CARRIER_IDS for participant in right)
        if not has_donor or not has_carrier:
            return False

        left_proteins = [
            participant for participant in left if participant.polymer_type == PolymerType.PROTEIN
        ]
        right_proteins = [
            participant for participant in right if participant.polymer_type == PolymerType.PROTEIN
        ]
        if len(left_proteins) != 1 or len(right_proteins) != 1:
            return False

        left_protein = left_proteins[0]
        right_protein = right_proteins[0]
        return left_protein.chebi_id != right_protein.chebi_id

    def _carrier_carbons(self, participant: Participant) -> int | None:
        """Return carrier carbon counts for released nucleotide products."""
        if participant.chebi_id in self.CARRIER_CARBONS:
            return self.CARRIER_CARBONS[participant.chebi_id]
        carbon_count = self._carbon_count(participant)
        if carbon_count and participant.chebi_id in self.CARRIER_CARBONS:
            return carbon_count
        return None

    def _ring_oxygen_count(self, participant: Participant) -> int:
        """Count ring oxygens as a proxy for cyclic glycosyl content."""
        mol = participant.get_mol()
        if mol is None or self.RING_OXYGEN_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(self.RING_OXYGEN_PATTERN))

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        """Count carbon atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        """Count nitrogen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        """Count phosphorus atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
