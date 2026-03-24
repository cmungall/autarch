"""Catalysis of the transfer of a pentosyl group from one compound (donor) to another (acceptor)."""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CO2,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_UDP,
)
from autarch.ontology.reaction import ReactionClass

CHEBI_DTDP = "CHEBI:16334"


class Pentosyltransferase(ReactionClass):
    """pentosyltransferase

    Catalysis of the transfer of a pentosyl group from one compound (donor) to another (acceptor).
    """

    GO_ID = "GO:0016763"
    EC_NUMBER_PREFIX = "2.4.2.-"
    SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS, CHEBI_CO2}
    CARRIER_CARBONS = {
        CHEBI_UDP: 9,
        CHEBI_GDP: 10,
        CHEBI_CDP: 9,
        CHEBI_DTDP: 10,
    }
    PENTOSE_RING_PATTERN = Chem.MolFromSmarts("[O;R]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for pentose transfer using carrier and carbon-balance constraints."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not pentosyl transfer",
            )

        chebi_ids = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if self._has_energy_coupling(chebi_ids):
            return ClassificationResult(
                is_member=False,
                explanation="ATP- or GTP-coupled chemistry indicates a different enzyme class",
            )
        if self._has_redox_pair(chebi_ids):
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P)-linked redox chemistry is not pentosyl transfer",
            )
        if CHEBI_O2 in chebi_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen chemistry is not pentosyl transfer",
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
        left = [participant for participant in left_participants if participant.chebi_id not in self.SPECTATOR_CHEBIS]
        right = [participant for participant in right_participants if participant.chebi_id not in self.SPECTATOR_CHEBIS]

        if self._is_nucleotide_pentose_transfer(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Pentosyltransferase: nucleotide carrier releases a pentose donor to a new acceptor",
            )

        if self._is_pentose_phosphorolysis(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Pentosyltransferase: phosphate or diphosphate accepts a transferred pentose group",
            )

        if self._is_phosphoribosyl_exchange(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Pentosyltransferase: a phosphorylated pentosyl donor exchanges the pentose between acceptors",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No pentosyl-transfer signature detected",
        )

    def _is_nucleotide_pentose_transfer(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Detect NDP- or dTDP-linked pentose transfer."""
        for carrier in right:
            carrier_carbons = self._carrier_carbons(carrier)
            if carrier_carbons is None:
                continue
            for donor in left:
                if not self._is_pentosylated_nucleotide_donor(donor, carrier_carbons):
                    continue
                for acceptor in left:
                    if acceptor is donor:
                        continue
                    if self._phosphorus_count(acceptor) > 0 or self._pentose_ring_count(acceptor) > 0:
                        continue
                    acceptor_carbons = self._carbon_count(acceptor)
                    if acceptor_carbons == 0:
                        continue
                    for product in right:
                        if product is carrier:
                            continue
                        if (
                            self._pentose_ring_count(product) >= 1
                            and self._carbon_count(product) == acceptor_carbons + 5
                        ):
                            return True
        return False

    def _is_pentose_phosphorolysis(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Detect nucleoside phosphorolysis and PRPP-forming exchange."""
        if not any(self._is_inorganic_phosphate(participant) for participant in left):
            return False

        if not any(self._is_free_pentose_phosphate(participant) for participant in right):
            return False

        for substrate in left:
            substrate_carbons = self._carbon_count(substrate)
            if self._pentose_ring_count(substrate) == 0 or substrate_carbons < 5:
                continue
            for leaving_group in right:
                if self._pentose_ring_count(leaving_group) > 0 or self._phosphorus_count(leaving_group) > 0:
                    continue
                if substrate_carbons == self._carbon_count(leaving_group) + 5:
                    return True
        return False

    def _is_phosphoribosyl_exchange(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Detect exchange of a phosphorylated pentose between acceptors."""
        for donor in left:
            if not self._is_phosphorylated_pentosylated_participant(donor):
                continue
            donor_carbons = self._carbon_count(donor)
            for leaving_group in right:
                if self._pentose_ring_count(leaving_group) > 0 or self._phosphorus_count(leaving_group) > 0:
                    continue
                if donor_carbons != self._carbon_count(leaving_group) + 5:
                    continue
                for acceptor in left:
                    if acceptor is donor:
                        continue
                    if self._phosphorus_count(acceptor) > 0 or self._pentose_ring_count(acceptor) > 0:
                        continue
                    acceptor_carbons = self._carbon_count(acceptor)
                    if acceptor_carbons == 0:
                        continue
                    for product in right:
                        if product is leaving_group:
                            continue
                        if (
                            self._is_phosphorylated_pentosylated_participant(product)
                            and self._carbon_count(product) == acceptor_carbons + 5
                        ):
                            return True
        return False

    def _is_pentosylated_nucleotide_donor(self, participant: Participant, carrier_carbons: int) -> bool:
        """Detect a nucleotide-linked donor carrying exactly one pentose more than the released carrier."""
        return (
            self._phosphorus_count(participant) >= 2
            and self._pentose_ring_count(participant) >= 1
            and self._nitrogen_count(participant) >= 1
            and self._carbon_count(participant) == carrier_carbons + 5
        )

    def _is_phosphorylated_pentosylated_participant(self, participant: Participant) -> bool:
        """Detect phosphoribosylated donors or products."""
        return (
            self._pentose_ring_count(participant) >= 1
            and self._phosphorus_count(participant) == 1
            and self._carbon_count(participant) >= 5
        )

    def _is_free_pentose_phosphate(self, participant: Participant) -> bool:
        """Detect a free pentose phosphate or diphosphate product."""
        return (
            self._pentose_ring_count(participant) >= 1
            and self._phosphorus_count(participant) >= 1
            and self._nitrogen_count(participant) == 0
            and self._carbon_count(participant) == 5
        )

    def _is_inorganic_phosphate(self, participant: Participant) -> bool:
        """Detect phosphate or diphosphate participants."""
        return self._carbon_count(participant) == 0 and self._phosphorus_count(participant) >= 1

    def _carrier_carbons(self, participant: Participant) -> int | None:
        """Get the carbon count for a released nucleotide carrier."""
        if participant.chebi_id in self.CARRIER_CARBONS:
            return self.CARRIER_CARBONS[participant.chebi_id]
        carbon_count = self._carbon_count(participant)
        if carbon_count and participant.chebi_id in {CHEBI_UDP, CHEBI_GDP, CHEBI_CDP, CHEBI_DTDP}:
            return carbon_count
        return None

    def _pentose_ring_count(self, participant: Participant) -> int:
        """Count ring oxygens as a proxy for cyclic pentose units."""
        mol = participant.get_mol()
        if mol is None or self.PENTOSE_RING_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(self.PENTOSE_RING_PATTERN))

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

    @staticmethod
    def _has_energy_coupling(chebi_ids: set[str]) -> bool:
        """Detect ATP- or GTP-coupled chemistry."""
        return (
            CHEBI_ATP in chebi_ids and bool({CHEBI_ADP, CHEBI_AMP} & chebi_ids)
        ) or (
            CHEBI_GTP in chebi_ids and CHEBI_GDP in chebi_ids
        )

    @staticmethod
    def _has_redox_pair(chebi_ids: set[str]) -> bool:
        """Detect NAD(P)-linked redox chemistry."""
        return (
            CHEBI_NAD_PLUS in chebi_ids and CHEBI_NADH in chebi_ids
        ) or (
            CHEBI_NADP_PLUS in chebi_ids and CHEBI_NADPH in chebi_ids
        )
