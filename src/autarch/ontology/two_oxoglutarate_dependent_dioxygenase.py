"""2-oxoglutarate-dependent dioxygenase.

Catalysis of the reaction: A + 2-oxoglutarate + O2 = B + succinate + CO2. This
is an oxidation-reduction (redox) reaction in which hydrogen or electrons are
transferred from 2-oxoglutarate and one other donor, and one atom of oxygen is
incorporated into each donor.
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


class TwoOxoglutarateDependentDioxygenase(Oxidoreductase):
    """2-oxoglutarate-dependent dioxygenase.

    Catalysis of the reaction: A + 2-oxoglutarate + O2 = B + succinate + CO2.
    This is an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from 2-oxoglutarate and one other donor, and one
    atom of oxygen is incorporated into each donor.
    """

    GO_ID = "GO:0016706"
    EC_NUMBER_PREFIX = "1.14.11.-"
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_FERREDOXIN_OXIDIZED,
        CHEBI_FERREDOXIN_REDUCED,
        CHEBI_H2O2,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_2_OXOGLUTARATE, CHEBI_O2, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_SUCCINATE, CHEBI_CO2, CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for the 2-oxoglutarate/succinate dioxygenase cochemistry."""
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        chebis = {
            participant.chebi_id
            for participant in left + right
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Nicotinamide, flavin, ferredoxin, or peroxide cofactors indicate a different oxygenase branch",
            )

        if not any(participant.chebi_id == CHEBI_2_OXOGLUTARATE for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="No 2-oxoglutarate cosubstrate in this reaction orientation",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen cosubstrate in this reaction orientation",
            )

        if not any(participant.chebi_id == CHEBI_SUCCINATE for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="No succinate coproduct detected",
            )

        if not any(participant.chebi_id == CHEBI_CO2 for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="No carbon dioxide coproduct detected",
            )

        left_reactive = [
            participant for participant in left if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant for participant in right if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive donor/product pair remains after removing 2-oxoglutarate cochemistry",
            )

        if not any(self._has_carbon(participant) for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No substantive organic donor remains after removing 2-oxoglutarate and dioxygen",
            )

        if not any(self._has_carbon(participant) for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No substantive organic product remains after removing succinate and carbon dioxide",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "2-oxoglutarate-dependent dioxygenase: dioxygen and 2-oxoglutarate "
                "are coupled to succinate and carbon dioxide formation"
            ),
        )

    @staticmethod
    def _has_carbon(participant: Participant) -> bool:
        """Check whether a participant has at least one carbon atom."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
