"""oxidoreductase acting on single donors with incorporation of molecular
oxygen incorporation of two atoms of oxygen.

Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
electrons are transferred from one donor, and two oxygen atoms is incorporated
into a donor.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_CO2,
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


class OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfTwoAtomsOfOxygen(
    Oxidoreductase
):
    """oxidoreductase acting on single donors with incorporation of molecular
    oxygen incorporation of two atoms of oxygen.

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from one donor, and two oxygen atoms is
    incorporated into a donor.
    """

    GO_ID = "GO:0016702"
    EC_NUMBER_PREFIX = "1.13.11.-"
    REDOX_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H_PLUS, CHEBI_H2O}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H_PLUS, CHEBI_H2O, CHEBI_H2O2, CHEBI_CO2}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for dioxygenase chemistry that incorporates both oxygen atoms."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        forward_result = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward_result.is_member:
            return forward_result

        reverse_result = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse_result.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse_result.explanation} (reverse reaction orientation)",
            )

        return forward_result

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen substrate in this reaction orientation",
            )

        if any(participant.chebi_id in self.REDOX_CHEBIS for participant in left_participants + right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P) chemistry indicates paired-donor oxygenation, not single-donor dioxygenation",
            )

        if any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Hydrogen peroxide product indicates oxidase chemistry, not incorporation of two oxygen atoms",
            )

        left_organic = [
            participant
            for participant in left_participants
            if self._is_carbon_containing(participant, self.LEFT_SPECTATOR_CHEBIS)
        ]
        if len(left_organic) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Single-donor dioxygenases should have exactly one substantive organic donor in one reaction orientation",
            )

        right_organic = [
            participant
            for participant in right_participants
            if self._is_carbon_containing(participant, self.RIGHT_SPECTATOR_CHEBIS)
        ]
        if not right_organic:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive organic oxygenated products detected",
            )

        oxygen_delta = self._oxygen_count(right_organic) - self._oxygen_count(left_organic)
        if oxygen_delta != 2:
            return ClassificationResult(
                is_member=False,
                explanation="Organic products do not gain two oxygen atoms",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Oxidoreductase acting on single donors with incorporation of molecular oxygen, "
                "incorporation of two atoms of oxygen"
            ),
        )

    @staticmethod
    def _is_carbon_containing(participant: Participant, excluded_chebis: set[str]) -> bool:
        """Check if a participant is a substantive carbon-containing species."""
        if participant.chebi_id in excluded_chebis:
            return False
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    @staticmethod
    def _oxygen_count(participants: list[Participant]) -> int:
        """Count oxygen atoms across participants."""
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            total += oxygen_count * max(participant.count, 1)
        return total
