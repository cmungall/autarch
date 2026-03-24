"""Catalysis of an oxidation-reduction (redox) reaction in which a sulfur-containing group acts as a hydrogen or electron donor and reduces disulfide.

This EC 1.8.4 branch captures sulfur-group redox reactions whose terminal
acceptor is a disulfide. Canonical examples include glutathione- and
thioredoxin-linked disulfide exchange, as well as sulfur oxygenations coupled to
reduction of a protein disulfide.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_FMN,
    CHEBI_FMNH2,
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

CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"


class OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on a sulfur group of donors, disulfide as acceptor

    Catalysis of an oxidation-reduction (redox) reaction in which a
    sulfur-containing group acts as a hydrogen or electron donor and reduces
    disulfide.

    The defining feature is reduction of a disulfide acceptor without external
    nicotinamide, flavin, peroxide, or dioxygen cochemistry. The partner sulfur
    donor is oxidized either to a new disulfide or to a more oxygenated sulfur
    state.
    """

    GO_ID = "GO:0016671"
    EC_NUMBER_PREFIX = "1.8.4.-"

    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_FERREDOXIN_OXIDIZED,
        CHEBI_FERREDOXIN_REDUCED,
    }
    SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    THIOL_PATTERN = Chem.MolFromSmarts("[S;X2H1,X1-;!$(S(=O));!$(SS)]")
    DISULFIDE_PATTERN = Chem.MolFromSmarts("[S;X2][S;X2]")
    SULFUR_OXYGEN_PATTERN = Chem.MolFromSmarts("[S;X3,X4](=[O;X1])")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for disulfide-acceptor sulfur redox chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sulfur-group redox chemistry",
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
                explanation="External nicotinamide, flavin, ferredoxin, peroxide, or dioxygen cofactors indicate another sulfur redox branch",
            )

        left_sulfur = self._sulfur_reactive(left_participants)
        right_sulfur = self._sulfur_reactive(right_participants)
        if len(left_sulfur) < 2 or len(right_sulfur) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Requires sulfur-bearing donor and acceptor chemistry on both sides of the reaction",
            )

        if not self._has_disulfide_acceptor_signature(left_sulfur, right_sulfur):
            return ClassificationResult(
                is_member=False,
                explanation="No disulfide acceptor is reduced in this reaction orientation",
            )

        if not self._has_partner_sulfur_oxidation(left_sulfur, right_sulfur):
            return ClassificationResult(
                is_member=False,
                explanation="No partner sulfur donor is oxidized as the disulfide acceptor is reduced",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sulfur-group oxidoreductase: a disulfide acceptor is reduced while a sulfur donor is oxidized",
        )

    def _sulfur_reactive(self, participants: list[Participant]) -> list[Participant]:
        return [
            participant
            for participant in participants
            if participant.chebi_id not in self.SPECTATOR_CHEBIS and self._sulfur_count(participant) > 0
        ]

    def _has_disulfide_acceptor_signature(
        self,
        left_sulfur: list[Participant],
        right_sulfur: list[Participant],
    ) -> bool:
        left_disulfides = sum(self._substructure_count(p, self.DISULFIDE_PATTERN) for p in left_sulfur)
        right_thiols = sum(self._substructure_count(p, self.THIOL_PATTERN) for p in right_sulfur)
        return left_disulfides >= 1 and right_thiols >= 2

    def _has_partner_sulfur_oxidation(
        self,
        left_sulfur: list[Participant],
        right_sulfur: list[Participant],
    ) -> bool:
        left_thiols = sum(self._substructure_count(p, self.THIOL_PATTERN) for p in left_sulfur)
        right_disulfides = sum(self._substructure_count(p, self.DISULFIDE_PATTERN) for p in right_sulfur)
        if left_thiols >= 2 and right_disulfides >= 1:
            return True

        left_s_oxygen = sum(self._substructure_count(p, self.SULFUR_OXYGEN_PATTERN) for p in left_sulfur)
        right_s_oxygen = sum(self._substructure_count(p, self.SULFUR_OXYGEN_PATTERN) for p in right_sulfur)
        return right_s_oxygen > left_s_oxygen

    @staticmethod
    def _substructure_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)

    @staticmethod
    def _sulfur_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16) * max(participant.count, 1)
