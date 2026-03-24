"""oxidoreductase acting on the CH-CH group of donors, iron-sulfur protein as
acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH group
acts as a hydrogen or electron donor and reduces an iron-sulfur protein.
"""

from __future__ import annotations

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_OXIDIZED_4FE4S = "CHEBI:33722"
CHEBI_REDUCED_4FE4S = "CHEBI:33723"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"


class OxidoreductaseActingOnTheCHCHGroupOfDonorsIronSulfurProteinAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase acting on the CH-CH group of donors, iron-sulfur protein as
    acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH
    group acts as a hydrogen or electron donor and reduces an iron-sulfur
    protein.
    """

    GO_ID = "GO:0016636"
    EC_NUMBER_PREFIX = "1.3.7.-"

    OXIDIZED_ACCEPTORS = {CHEBI_OXIDIZED_4FE4S, CHEBI_FERREDOXIN_OXIDIZED}
    REDUCED_ACCEPTORS = {CHEBI_REDUCED_4FE4S, CHEBI_FERREDOXIN_REDUCED}
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        "CHEBI:58210",  # FMN
        "CHEBI:58307",  # FMNH2
        CHEBI_O2,
        CHEBI_H2O2,
    }
    SPECTATOR_CHEBIS = OXIDIZED_ACCEPTORS | REDUCED_ACCEPTORS | {CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for CH-CH redox with iron-sulfur acceptors."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH-CH iron-sulfur redox chemistry",
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
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide, flavin, or oxygen cofactors indicate another CH-CH oxidoreductase branch",
            )

        left_acceptor_count = self._count_chebis(left_participants, self.OXIDIZED_ACCEPTORS)
        right_acceptor_count = self._count_chebis(right_participants, self.REDUCED_ACCEPTORS)
        if left_acceptor_count == 0 or left_acceptor_count != right_acceptor_count:
            return ClassificationResult(
                is_member=False,
                explanation="Missing oxidized-to-reduced iron-sulfur acceptor conversion",
            )

        left_reactive = self._reactive_carbon_participants(left_participants)
        right_reactive = self._reactive_carbon_participants(right_participants)
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive carbon-containing substrate/product pair remained after removing iron-sulfur acceptors and solvent spectators",
            )

        if not self._has_ch_ch_redox_pair(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No CH-CH substrate/product pair showed the expected unsaturation change",
            )

        return ClassificationResult(
            is_member=True,
            explanation="CH-CH oxidoreductase: an oxidized iron-sulfur acceptor is reduced while a CH-CH substrate changes unsaturation",
        )

    @staticmethod
    def _count_chebis(participants: list[Participant], chebis: set[str]) -> int:
        return sum(
            max(participant.count, 1)
            for participant in participants
            if participant.chebi_id in chebis
        )

    def _reactive_carbon_participants(
        self,
        participants: list[Participant],
    ) -> list[Participant]:
        return [
            participant
            for participant in participants
            if participant.chebi_id not in self.SPECTATOR_CHEBIS
            and participant.get_mol() is not None
            and self._carbon_count(participant) > 0
        ]

    def _has_ch_ch_redox_pair(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        for substrate in left_reactive:
            substrate_counts = self._heavy_atom_counts(substrate)
            substrate_cc_double = self._cc_double_bond_count(substrate)
            for product in right_reactive:
                if substrate_counts != self._heavy_atom_counts(product):
                    continue
                if substrate_cc_double != self._cc_double_bond_count(product):
                    return True
        return False

    @staticmethod
    def _heavy_atom_counts(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() != 1
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _cc_double_bond_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        count = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                count += 1
        return count
