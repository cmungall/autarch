"""Other lyase reaction classification.

EC 4.99: Other lyases that do not fit into the standard lyase subcategories.
These include enzymes like ferrochelatase that catalyze unusual bond
cleavage/formation reactions.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase


class OtherLyase(Lyase):
    """other lyase

    EC 4.99.x.x includes:
    - Ferrochelatase (4.99.1.1): inserts Fe2+ into protoporphyrin IX
    - Sirohydrochlorin cobaltochelatase
    - Other chelatases

    Examples:
    - Protoporphyrin IX + Fe2+ = protoheme + 2 H+
    """

    GO_ID = "GO:0016829"  # lyase activity (same as parent - catch-all)
    EC_NUMBER_PREFIX = "4.99.-.-"

    OTHER_LYASE_PATTERNS = [
        "ferrochelatase", "chelatase",
        "protoporphyrin", "cobaltochelatase",
        "sirohydrochlorin",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an 'other' lyase (EC 4.99).

        Strategy:
        1. Check for specific EC 4.99 patterns (chelatases)
        2. Otherwise delegate to parent Lyase
        """
        label_lower = reaction.label.lower() if reaction.label else ""

        has_other_pattern = any(
            pattern in label_lower
            for pattern in self.OTHER_LYASE_PATTERNS
        )

        if has_other_pattern:
            return ClassificationResult(
                is_member=True,
                explanation="Other lyase: chelatase/metal insertion"
            )

        # Delegate to parent for general lyase check
        parent_result = super().check_membership_impl(reaction)
        return parent_result
