"""Pectate lyase reaction classification."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O
from autarch.ontology.lyase import Lyase


class PectateLyase(Lyase):
    """pectate lyase"""

    GO_ID = "GO:0030570"  # pectate lyase activity
    EC_NUMBER_PREFIX = "4.2.2.2"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for eliminative cleavage of pectate (not hydrolysis)."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a lyase: {parent_result.explanation}",
            )

        has_water_reactant = any(
            p.chebi_id == CHEBI_H2O for p in reaction.left_participants
        )
        if has_water_reactant:
            return ClassificationResult(
                is_member=False,
                explanation="Water reactant indicates hydrolysis, not pectate lyase",
            )

        left_text = " ".join((p.name or "") for p in reaction.left_participants).lower()
        right_text = " ".join((p.name or "") for p in reaction.right_participants).lower()
        label_text = (reaction.label or "").lower()

        has_pectate_substrate = any(
            keyword in left_text or keyword in label_text
            for keyword in ("pectate", "galacturonan", "polygalacturon")
        )
        if not has_pectate_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No pectate/galacturonan substrate signature",
            )

        has_unsaturated_uronate_product = any(
            keyword in right_text or keyword in label_text
            for keyword in ("4-deoxy", "enuronosyl", "5-dehydro-4-deoxy")
        )
        if not has_unsaturated_uronate_product:
            return ClassificationResult(
                is_member=False,
                explanation="No unsaturated uronate product signature",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Pectate lyase: eliminative cleavage to unsaturated uronates",
        )
