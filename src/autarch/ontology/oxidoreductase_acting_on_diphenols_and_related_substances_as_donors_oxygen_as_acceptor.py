"""oxidoreductase acting on diphenols and related substances as donors, oxygen as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a diphenol, or
related compound, acts as a hydrogen or electron donor and reduces oxygen.
"""

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


class OxidoreductaseActingOnDiphenolsAndRelatedSubstancesAsDonorsOxygenAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase acting on diphenols and related substances as donors, oxygen as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a diphenol,
    or related compound, acts as a hydrogen or electron donor and reduces
    oxygen.
    """

    GO_ID = "GO:0016682"
    EC_NUMBER_PREFIX = "1.10.3.-"
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H2O2, CHEBI_H_PLUS}
    AROMATIC_HYDROXY_PATTERN = Chem.MolFromSmarts("[c][OX2H,OX1-]")
    ENEDIOL_PATTERN = Chem.MolFromSmarts("[C]=[C]([O;H1,-])[O;H1,-]")
    ENEDIOL_PATTERN_ALT = Chem.MolFromSmarts("[C]([O;H1,-])=[C]([O;H1,-])")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")
    RELATED_DONOR_CHEBIS = {"CHEBI:36559", "CHEBI:38290"}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for oxygen-acceptor oxidation of diphenols or related donors."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not diphenol oxidation",
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
        if any(
            participant.chebi_id in self.EXCLUDED_CHEBIS
            for participant in left_participants + right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide or flavin cofactors indicate a different diphenol oxidoreductase branch",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen acceptor in this reaction orientation",
            )

        if not any(
            participant.chebi_id in {CHEBI_H2O, CHEBI_H2O2}
            for participant in right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen-acceptor diphenol oxidation should yield water or hydrogen peroxide",
            )

        left_organics = self._organic_participants(left_participants, self.LEFT_SPECTATOR_CHEBIS)
        right_organics = self._organic_participants(right_participants, self.RIGHT_SPECTATOR_CHEBIS)
        donor_candidates = [
            participant
            for participant in left_organics
            if self._diphenol_score(participant) >= 2
            or participant.chebi_id in self.RELATED_DONOR_CHEBIS
        ]
        if not donor_candidates or not right_organics:
            return ClassificationResult(
                is_member=False,
                explanation="No diphenol or related donor substrate detected",
            )

        if self._has_ascorbate_cossubstrate_oxygenation(left_organics, donor_candidates, right_organics):
            return ClassificationResult(
                is_member=False,
                explanation="Ascorbate-like donor is acting as a cosubstrate while another organic substrate undergoes oxygenation",
            )

        if not self._has_oxidized_pair(donor_candidates, right_organics):
            return ClassificationResult(
                is_member=False,
                explanation="Organic products do not show oxidation of a diphenol-like donor scaffold",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Diphenol oxidoreductase with oxygen as acceptor",
        )

    @staticmethod
    def _organic_participants(
        participants: list[Participant],
        excluded_chebis: set[str],
    ) -> list[Participant]:
        organics: list[Participant] = []
        for participant in participants:
            if participant.chebi_id in excluded_chebis:
                continue
            mol = participant.get_mol()
            if mol is None:
                continue
            if any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
                organics.append(participant)
        return organics

    def _diphenol_score(self, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        aromatic_hydroxyls = len(mol.GetSubstructMatches(self.AROMATIC_HYDROXY_PATTERN))
        enediols = len(mol.GetSubstructMatches(self.ENEDIOL_PATTERN)) + len(
            mol.GetSubstructMatches(self.ENEDIOL_PATTERN_ALT)
        )
        return aromatic_hydroxyls + (2 * enediols)

    def _diphenol_total(self, participants: list[Participant]) -> int:
        return sum(self._diphenol_score(participant) * max(participant.count, 1) for participant in participants)

    def _carbonyl_count(self, participants: list[Participant]) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += len(mol.GetSubstructMatches(self.CARBONYL_PATTERN)) * max(participant.count, 1)
        return total

    def _has_oxidized_pair(
        self,
        donors: list[Participant],
        products: list[Participant],
    ) -> bool:
        for donor in donors:
            for product in products:
                if self._carbon_count(donor) != self._carbon_count(product):
                    continue
                if self._ring_count(donor) != self._ring_count(product):
                    continue
                if self._carbonyl_count([product]) > self._carbonyl_count([donor]):
                    return True
                if self._diphenol_score(product) < self._diphenol_score(donor):
                    return True
                if donor.chebi_id == "CHEBI:38290" and product.chebi_id == "CHEBI:59513":
                    return True
        return False

    def _has_ascorbate_cossubstrate_oxygenation(
        self,
        left_participants: list[Participant],
        donor_candidates: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        if not any(participant.chebi_id == "CHEBI:38290" for participant in donor_candidates):
            return False
        for left in left_participants:
            if left.chebi_id == "CHEBI:38290":
                continue
            for right in right_participants:
                if self._carbon_count(left) != self._carbon_count(right):
                    continue
                if right.is_same_molecule(left):
                    continue
                left_oxygen = self._oxygen_count(left)
                right_oxygen = self._oxygen_count(right)
                if right_oxygen > left_oxygen:
                    return True
        return False

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @staticmethod
    def _ring_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return mol.GetRingInfo().NumRings()
