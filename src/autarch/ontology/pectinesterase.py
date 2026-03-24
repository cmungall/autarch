"""Pectinesterase reaction classification."""

from autarch.datamodel import ClassificationResult, Reaction, PolymerType
from autarch.ontology.reaction import ReactionClass


CHEBI_WATER = "CHEBI:15377"
CHEBI_H_PLUS = "CHEBI:15378"
CHEBI_METHANOL = "CHEBI:17790"
CHEBI_GALACTURONOSYL_METHYL_ESTER = "CHEBI:140522"
CHEBI_GALACTURONOSYL = "CHEBI:140523"


class Pectinesterase(ReactionClass):
    """pectinesterase"""

    GO_ID = "GO:0030599"  # pectinesterase activity
    EC_NUMBERS = ["3.1.1.11"]
    EC_NUMBER_PREFIX = None

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for galacturonan methyl-ester hydrolysis."""
        left_polymers = self._find_galacturonan_polymers(reaction.left_participants)
        right_polymers = self._find_galacturonan_polymers(reaction.right_participants)

        if not left_polymers or not right_polymers:
            return ClassificationResult(
                is_member=False,
                explanation="No galacturonan polymer found on both sides",
            )

        left_has_ester = any(self._is_ester(p) for p in left_polymers)
        right_has_ester = any(self._is_ester(p) for p in right_polymers)
        left_has_deesterified = any(self._is_deesterified(p) for p in left_polymers)
        right_has_deesterified = any(self._is_deesterified(p) for p in right_polymers)

        forward = left_has_ester and right_has_deesterified
        reverse = right_has_ester and left_has_deesterified
        if not (forward or reverse):
            return ClassificationResult(
                is_member=False,
                explanation="No methyl-ester to pectate conversion pattern",
            )

        methanol_side = reaction.right_participants if forward else reaction.left_participants
        if not self._has_chebi(methanol_side, CHEBI_METHANOL):
            return ClassificationResult(
                is_member=False,
                explanation="No methanol on ester-cleavage side",
            )

        water_side = reaction.left_participants if forward else reaction.right_participants
        if not self._has_chebi(water_side, CHEBI_WATER):
            return ClassificationResult(
                is_member=False,
                explanation="No water on hydrolysis side",
            )

        if forward:
            return ClassificationResult(
                is_member=True,
                explanation="Pectinesterase: pectin methyl-ester hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Pectinesterase (reverse): pectate esterification",
        )

    def _find_galacturonan_polymers(self, participants):
        """Find valid galacturonan polymer participants."""
        GALACTURONAN_CHEBI = {
            "CHEBI:28716",
            CHEBI_GALACTURONOSYL_METHYL_ESTER,
            CHEBI_GALACTURONOSYL,
        }
        polymers = []
        for p in participants:
            if p.chebi_id in {CHEBI_WATER, CHEBI_H_PLUS}:
                continue
            if p.polymer_type == PolymerType.GALACTURONAN or p.chebi_id in GALACTURONAN_CHEBI:
                polymers.append(p)
        return polymers

    def _is_ester(self, participant) -> bool:
        """Check if participant is the methyl-ester polymer form."""
        monomer = (participant.monomer or "").lower()
        name = (participant.name or "").lower()
        return (
            participant.chebi_id == CHEBI_GALACTURONOSYL_METHYL_ESTER
            or "methyl ester" in monomer
            or "methyl ester" in name
        )

    def _is_deesterified(self, participant) -> bool:
        """Check if participant is the de-esterified polymer form."""
        if self._is_ester(participant):
            return False
        monomer = (participant.monomer or "").lower()
        name = (participant.name or "").lower()
        return (
            participant.chebi_id == CHEBI_GALACTURONOSYL
            or "galacturonosyl" in monomer
            or "galacturonosyl" in name
        )

    def _has_chebi(self, participants, chebi_id: str) -> bool:
        return any(p.chebi_id == chebi_id for p in participants)
