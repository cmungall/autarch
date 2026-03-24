"""oxidoreductase activity, acting on paired donors, with incorporation or
reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of two
atoms of oxygen into one donor.

Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
electrons are transferred from NADH or NADPH and one other donor, and two atoms
of oxygen are incorporated into one donor.
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
    CHEBI_NH4,
    CHEBI_NH3,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"
CHEBI_IRON_II = "CHEBI:29033"
CHEBI_IRON_III = "CHEBI:29034"


class OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfTwoAtomsOfOxygenIntoOneDonor(
    Oxidoreductase
):
    """oxidoreductase activity, acting on paired donors, with incorporation or
    reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of
    two atoms of oxygen into one donor.

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from NADH or NADPH and one other donor, and two
    atoms of oxygen are incorporated into one donor.
    """

    GO_ID = "GO:0016708"
    EC_NUMBER_PREFIX = "1.14.12.-"

    NAD_REDUCED_CHEBIS = {CHEBI_NADH, CHEBI_NADPH}
    NAD_OXIDIZED_CHEBIS = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    EXCLUDED_CHEBIS = {
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_H2O2,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH}
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_H_PLUS,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_CO2,
        CHEBI_NH3,
        CHEBI_NH4,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for NAD(P)H-linked dioxygenation chemistry."""
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
                explanation="Flavin or peroxide chemistry indicates a different oxygenase branch",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen substrate in this reaction orientation",
            )

        if any(participant.chebi_id == CHEBI_H2O for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Water coproduct indicates one-atom oxygenation rather than dioxygenation",
            )

        if not self._has_redox_pair(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing NAD(P)H-to-NAD(P)+ conversion",
            )

        if self._has_2_oxoglutarate_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate/succinate cochemistry indicates a different oxygenase branch",
            )

        if self._has_ferredoxin_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Reduced iron-sulfur protein chemistry indicates a different oxygenase branch",
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
        if len(left_organic) != 1 or not right_organic:
            return ClassificationResult(
                is_member=False,
                explanation="Dioxygenation requires one substantive organic donor and at least one substantive organic product",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Oxidoreductase acting on paired donors with incorporation or reduction "
                "of molecular oxygen, NAD(P)H as one donor, and incorporation of two atoms "
                "of oxygen into one donor"
            ),
        )

    @classmethod
    def _has_redox_pair(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for reduced-to-oxidized nicotinamide conversion."""
        has_reduced = any(
            participant.chebi_id in cls.NAD_REDUCED_CHEBIS for participant in left_participants
        )
        has_oxidized = any(
            participant.chebi_id in cls.NAD_OXIDIZED_CHEBIS for participant in right_participants
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
    def _has_ferredoxin_branch(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for reduced iron-sulfur donor chemistry."""
        reduced_chebis = {CHEBI_FERREDOXIN_REDUCED, CHEBI_IRON_II}
        oxidized_chebis = {CHEBI_FERREDOXIN_OXIDIZED, CHEBI_IRON_III}
        has_reduced = any(participant.chebi_id in reduced_chebis for participant in left_participants)
        has_oxidized = any(participant.chebi_id in oxidized_chebis for participant in right_participants)
        return has_reduced and has_oxidized

    @staticmethod
    def _has_carbon(participant: Participant) -> bool:
        """Check whether a participant is a substantive carbon-containing species."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
