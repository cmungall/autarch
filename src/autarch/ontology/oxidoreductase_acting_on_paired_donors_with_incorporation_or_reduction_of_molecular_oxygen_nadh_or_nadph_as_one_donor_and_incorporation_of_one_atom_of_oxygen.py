"""oxidoreductase activity, acting on paired donors, with incorporation or
reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one
atom of oxygen.

Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
electrons are transferred from NADH or NADPH and one other donor, and one atom
of oxygen is incorporated into one donor.
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
CHEBI_SUCCINATE = "CHEBI:15741"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"


class OxidoreductaseActingOnPairedDonorsWithIncorporationOrReductionOfMolecularOxygenNADHOrNADPHAsOneDonorAndIncorporationOfOneAtomOfOxygen(
    Oxidoreductase
):
    """oxidoreductase activity, acting on paired donors, with incorporation or
    reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of
    one atom of oxygen.

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from NADH or NADPH and one other donor, and one
    atom of oxygen is incorporated into one donor.
    """

    GO_ID = "GO:0016709"
    EC_NUMBER_PREFIX = "1.14.13.-"
    NAD_REDUCED_CHEBIS = {CHEBI_NADH, CHEBI_NADPH}
    NAD_OXIDIZED_CHEBIS = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    FLAVIN_REDUCED_CHEBIS = {CHEBI_FADH2, CHEBI_FMNH2}
    FLAVIN_OXIDIZED_CHEBIS = {CHEBI_FAD, CHEBI_FMN}
    LEFT_SPECTATOR_CHEBIS = {
        CHEBI_O2,
        CHEBI_H_PLUS,
        CHEBI_H2O,
        CHEBI_NADH,
        CHEBI_NADPH,
        CHEBI_FADH2,
        CHEBI_FMNH2,
    }
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H2O2,
        CHEBI_H_PLUS,
        CHEBI_CO2,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_FAD,
        CHEBI_FMN,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for paired-donor oxygenase chemistry within the GO:0016709 branch."""
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

        if any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Hydrogen peroxide product indicates oxidase chemistry, not paired-donor oxygenation",
            )

        if self._has_ferredoxin_pair(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Ferredoxin/adrenodoxin chemistry indicates a different monooxygenase branch",
            )

        if self._has_2_oxoglutarate_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate/succinate cochemistry indicates a different oxygenase branch",
            )

        has_nad_pair = self._has_redox_pair(
            left_participants,
            right_participants,
            self.NAD_REDUCED_CHEBIS,
            self.NAD_OXIDIZED_CHEBIS,
        )
        has_flavin_pair = self._has_redox_pair(
            left_participants,
            right_participants,
            self.FLAVIN_REDUCED_CHEBIS,
            self.FLAVIN_OXIDIZED_CHEBIS,
        )

        if not (has_nad_pair or has_flavin_pair):
            return ClassificationResult(
                is_member=False,
                explanation="Missing paired-donor redox cofactor conversion",
            )

        if has_nad_pair and has_flavin_pair:
            return ClassificationResult(
                is_member=False,
                explanation="Mixed NAD(P)H/flavin relay chemistry indicates a different oxygenase branch",
            )

        has_water_product = any(
            participant.chebi_id == CHEBI_H2O for participant in right_participants
        )
        has_nadh = any(participant.chebi_id == CHEBI_NADH for participant in left_participants)
        has_nadph = any(participant.chebi_id == CHEBI_NADPH for participant in left_participants)
        if has_nad_pair and has_nadh and not has_nadph and not has_water_product:
            return ClassificationResult(
                is_member=False,
                explanation="NADH-dependent oxygenation without water product is more consistent with paired-donor dioxygenase chemistry",
            )

        substrate = self._dominant_organic_participant(
            left_participants,
            self.LEFT_SPECTATOR_CHEBIS,
        )
        product = self._dominant_organic_participant(
            right_participants,
            self.RIGHT_SPECTATOR_CHEBIS,
        )
        if substrate is None or product is None:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive donor/product pair detected",
            )

        if self._oxygen_count(product) < self._oxygen_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Dominant organic product loses oxygen rather than retaining/incorporating it",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Oxidoreductase acting on paired donors with incorporation or reduction "
                "of molecular oxygen and incorporation of one atom of oxygen"
            ),
        )

    @classmethod
    def _has_redox_pair(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
        reduced_chebis: set[str],
        oxidized_chebis: set[str],
    ) -> bool:
        """Check for reduced-to-oxidized cofactor conversion across the reaction."""
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
        """Check for ferredoxin/adrenodoxin-mediated oxygenase chemistry."""
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
        """Check for 2-oxoglutarate/succinate oxygenase chemistry."""
        has_2_oxoglutarate = any(
            participant.chebi_id == CHEBI_2_OXOGLUTARATE for participant in left_participants
        )
        has_succinate = any(
            participant.chebi_id == CHEBI_SUCCINATE for participant in right_participants
        )
        return has_2_oxoglutarate and has_succinate

    @staticmethod
    def _dominant_organic_participant(
        participants: list[Participant],
        excluded_chebis: set[str],
    ) -> Participant | None:
        """Return the largest substantive carbon-containing participant."""
        best_participant: Participant | None = None
        best_signature = (-1, -1)
        for participant in participants:
            if participant.chebi_id in excluded_chebis:
                continue
            mol = participant.get_mol()
            if mol is None:
                continue
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            if carbon_count == 0:
                continue
            oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            signature = (carbon_count, oxygen_count)
            if signature > best_signature:
                best_participant = participant
                best_signature = signature
        return best_participant

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        """Count oxygen atoms in a participant."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
