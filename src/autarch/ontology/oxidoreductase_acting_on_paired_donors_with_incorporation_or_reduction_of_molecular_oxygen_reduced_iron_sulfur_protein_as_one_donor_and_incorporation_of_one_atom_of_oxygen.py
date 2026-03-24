"""oxidoreductase acting on paired donors with incorporation or reduction of molecular oxygen reduced iron-sulfur protein as one donor and incorporation of one atom of oxygen.

Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
electrons are transferred from reduced iron-sulfur protein and one other donor,
and one atom of oxygen is incorporated into one donor.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_CO2,
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
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"
CHEBI_HIGH_NUCLEARITY_IRON_SULFUR_OXIDIZED = "CHEBI:33722"
CHEBI_HIGH_NUCLEARITY_IRON_SULFUR_REDUCED = "CHEBI:33723"
CHEBI_IRON_II = "CHEBI:29033"
CHEBI_IRON_III = "CHEBI:29034"


class OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenReducedIronSulfurProteinAsOneDonorAndIncorporationOfOneAtomOfOxygen(
    Oxidoreductase
):
    """oxidoreductase acting on paired donors with incorporation or reduction of molecular oxygen reduced iron-sulfur protein as one donor and incorporation of one atom of oxygen.

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from reduced iron-sulfur protein and one other
    donor, and one atom of oxygen is incorporated into one donor.
    """

    GO_ID = "GO:0016713"
    EC_NUMBER_PREFIX = "1.14.15.-"

    REDUCED_IRON_SULFUR_CHEBIS = {
        CHEBI_FERREDOXIN_REDUCED,
        CHEBI_HIGH_NUCLEARITY_IRON_SULFUR_REDUCED,
        CHEBI_IRON_II,
    }
    OXIDIZED_IRON_SULFUR_CHEBIS = {
        CHEBI_FERREDOXIN_OXIDIZED,
        CHEBI_HIGH_NUCLEARITY_IRON_SULFUR_OXIDIZED,
        CHEBI_IRON_III,
    }
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_H2O2,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H_PLUS} | REDUCED_IRON_SULFUR_CHEBIS
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS, CHEBI_CO2} | OXIDIZED_IRON_SULFUR_CHEBIS

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for reduced iron-sulfur protein-dependent monooxygenation."""
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
        """Evaluate one reaction orientation."""
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Nicotinamide, flavin, or peroxide chemistry indicates a different oxygenase branch",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen substrate in this reaction orientation",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="One-atom oxygenation branch requires water as the second oxygen sink",
            )

        if not self._has_iron_sulfur_pair(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing reduced-to-oxidized iron-sulfur donor conversion",
            )

        if self._has_2_oxoglutarate_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate/succinate cochemistry indicates a different oxygenase branch",
            )

        left_organic = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS and self._has_carbon(participant)
        ]
        right_organic = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS and self._has_carbon(participant)
        ]
        if not left_organic or not right_organic:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive donor/product pair detected",
            )

        if self._oxygen_count(right_organic) <= self._oxygen_count(left_organic):
            return ClassificationResult(
                is_member=False,
                explanation="Organic products do not gain oxygen overall, which is inconsistent with monooxygenation",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Oxidoreductase acting on paired donors with incorporation or reduction "
                "of molecular oxygen, reduced iron-sulfur protein as one donor, and "
                "incorporation of one atom of oxygen"
            ),
        )

    @classmethod
    def _has_iron_sulfur_pair(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for reduced-to-oxidized iron-sulfur donor conversion."""
        has_reduced = any(
            participant.chebi_id in cls.REDUCED_IRON_SULFUR_CHEBIS for participant in left_participants
        )
        has_oxidized = any(
            participant.chebi_id in cls.OXIDIZED_IRON_SULFUR_CHEBIS for participant in right_participants
        )
        return has_reduced and has_oxidized

    @staticmethod
    def _has_2_oxoglutarate_branch(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for 2-oxoglutarate-dependent oxygenase cochemistry."""
        has_2_oxoglutarate = any(
            participant.chebi_id == CHEBI_2_OXOGLUTARATE for participant in left_participants
        )
        has_succinate = any(
            participant.chebi_id == CHEBI_SUCCINATE for participant in right_participants
        )
        return has_2_oxoglutarate and has_succinate

    @staticmethod
    @staticmethod
    def _has_carbon(participant: Participant) -> bool:
        """Check whether a participant has at least one carbon atom."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    @staticmethod
    def _oxygen_count(participants: list[Participant]) -> int:
        """Count oxygen atoms across a participant list."""
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            total += oxygen_count * max(participant.count, 1)
        return total
