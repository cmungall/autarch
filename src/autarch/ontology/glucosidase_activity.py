"""glucosidase activity.

Catalysis of hydrolysis of a glucosidic bond with release of glucose or
glucose-phosphate products.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase_hydrolyzing_o_glycosyl_compounds import (
    HydrolaseHydrolyzingOGlycosylCompounds,
)
from autarch.ontology.lipid_utils import is_glucose_like


class GlucosidaseActivity(HydrolaseHydrolyzingOGlycosylCompounds):
    """glucosidase activity.

    Catalysis of hydrolysis of a glucosidic bond with release of glucose or
    glucose-phosphate products.
    """

    GO_ID = "GO:0015926"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        if not any(is_glucose_like(participant) for participant in reaction.right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No glucose-like product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Glucosidase activity: hydrolysis of a glucosidic bond with glucose release",
        )
