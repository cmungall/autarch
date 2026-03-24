"""oxidoreductase acting on paired donors with incorporation or reduction of molecular oxygen reduced flavin or flavoprotein as one donor and incorporation of one atom of oxygen.

Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or electrons are transferred from reduced flavin or flavoprotein and one other donor, and one atom of oxygen is incorporated into one donor.
"""

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

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:15741"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"


class OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedFlavinOrFlavoproteinAsOneDonorAndIncorporationOfOneAtomOfOxygen(
    Oxidoreductase
):
    """oxidoreductase acting on paired donors with incorporation or reduction of molecular oxygen reduced flavin or flavoprotein as one donor and incorporation of one atom of oxygen.

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from reduced flavin or flavoprotein and one other
    donor, and one atom of oxygen is incorporated into one donor.
    """

    GO_ID = "GO:0016712"
    EC_NUMBER_PREFIX = "1.14.14.-"

    REDUCED_FLAVIN_CHEBIS = {CHEBI_FADH2, CHEBI_FMNH2}
    OXIDIZED_FLAVIN_CHEBIS = {CHEBI_FAD, CHEBI_FMN}
    LEFT_SPECTATOR_CHEBIS = {
        CHEBI_O2,
        CHEBI_H_PLUS,
        CHEBI_FADH2,
        CHEBI_FMNH2,
        CHEBI_NADH,
        CHEBI_NADPH,
    }
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H2O2,
        CHEBI_H_PLUS,
        CHEBI_FAD,
        CHEBI_FMN,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for flavin-dependent monooxygenation."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
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
        """Evaluate one reaction direction."""
        o2_count = self._total_count(left_participants, CHEBI_O2)
        if o2_count == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen substrate in this reaction orientation",
            )

        if self._total_count(right_participants, CHEBI_H2O2) > 0:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrogen peroxide product indicates oxidase chemistry",
            )

        if self._total_count(right_participants, CHEBI_H2O) < o2_count:
            return ClassificationResult(
                is_member=False,
                explanation="Monooxygenation branch requires at least one water product per dioxygen",
            )

        if not self._has_redox_pair(
            left_participants,
            right_participants,
            self.REDUCED_FLAVIN_CHEBIS,
            self.OXIDIZED_FLAVIN_CHEBIS,
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Missing reduced-flavin to oxidized-flavin conversion",
            )

        if self._has_ferredoxin_pair(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Ferredoxin chemistry indicates a different oxygenase branch",
            )

        if self._has_2_oxoglutarate_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate/succinate cochemistry indicates a different oxygenase branch",
            )

        left_organics = self._organic_participants(
            left_participants,
            self.LEFT_SPECTATOR_CHEBIS,
        )
        right_organics = self._organic_participants(
            right_participants,
            self.RIGHT_SPECTATOR_CHEBIS,
        )
        if not left_organics or not right_organics:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive donor/product pair detected",
            )

        if self._total_oxygen(right_organics) <= self._total_oxygen(left_organics):
            return ClassificationResult(
                is_member=False,
                explanation="Organic products do not gain oxygen overall",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Oxidoreductase acting on paired donors with reduced flavin or "
                "flavoprotein as one donor and incorporation of one atom of oxygen"
            ),
        )

    @staticmethod
    def _total_count(participants: list[Participant], chebi_id: str) -> int:
        return sum(participant.count for participant in participants if participant.chebi_id == chebi_id)

    @staticmethod
    def _has_redox_pair(
        left_participants: list[Participant],
        right_participants: list[Participant],
        reduced_chebis: set[str],
        oxidized_chebis: set[str],
    ) -> bool:
        has_reduced = any(participant.chebi_id in reduced_chebis for participant in left_participants)
        has_oxidized = any(
            participant.chebi_id in oxidized_chebis for participant in right_participants
        )
        return has_reduced and has_oxidized

    @staticmethod
    def _has_ferredoxin_pair(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        has_reduced = any(
            participant.chebi_id == CHEBI_FERREDOXIN_REDUCED for participant in left_participants
        )
        has_oxidized = any(
            participant.chebi_id == CHEBI_FERREDOXIN_OXIDIZED for participant in right_participants
        )
        return has_reduced and has_oxidized

    @staticmethod
    def _has_2_oxoglutarate_branch(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        has_2_oxoglutarate = any(
            participant.chebi_id == CHEBI_2_OXOGLUTARATE for participant in left_participants
        )
        has_succinate = any(
            participant.chebi_id == CHEBI_SUCCINATE for participant in right_participants
        )
        return has_2_oxoglutarate and has_succinate

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

    @staticmethod
    def _total_oxygen(participants: list[Participant]) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            total += oxygen_atoms * max(participant.count, 1)
        return total
