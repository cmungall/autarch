"""beta-glucosidase activity."""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.glucosidase_activity import GlucosidaseActivity


class BetaGlucosidaseActivity(GlucosidaseActivity):
    """beta-glucosidase activity."""

    GO_ID = "GO:0008422"
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = None

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent
        return ClassificationResult(is_member=True, explanation="Beta-glucosidase activity: glucoside hydrolysis with glucose release")
