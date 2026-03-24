"""Intramolecular lyase (cyclase) reaction classification.

EC 5.5: Intramolecular lyases catalyze reactions where a bond is broken
within a molecule to form a ring (cyclization). Despite the name "lyase",
these are classified as isomerases (EC 5) because they rearrange a single
molecule without net bond cleavage.

Pattern: linear substrate → cyclic product (1→1 transformation)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase


class IntramolecularLyase(Isomerase):
    """intramolecular lyase

    EC 5.5.x.x includes:
    - Muconate cycloisomerase (5.5.1.1): muconate → muconolactone
    - Chorismate mutase (5.5.1.5) - not a true mutase
    - Strictosidine synthase-like cyclases
    - Various cyclization reactions forming lactones, lactams, or carbocycles

    Examples:
    - cis,cis-Muconate = muconolactone
    - Chorismate = prephenate
    """

    GO_ID = "GO:0016872"  # intramolecular lyase activity
    EC_NUMBER_PREFIX = "5.5.-.-"

    CYCLASE_LABEL_PATTERNS = [
        "cycloisomerase", "cyclase",
        "muconate", "muconolactone",
        "chorismate", "prephenate",
        "isopentenyl", "dimethylallyl",
        "copalyl", "ent-kaurene",
        "squalene-hopene",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an intramolecular lyase (cyclase).

        Strategy:
        1. Must be 1→1 stoichiometry (isomerization/cyclization)
        2. Look for cyclase/cycloisomerase label patterns
        """
        # Must be 1→1 (intramolecular rearrangement forming a ring)
        if len(reaction.left_participants) != 1 or len(reaction.right_participants) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Not 1→1 stoichiometry - not intramolecular lyase"
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        has_cyclase_pattern = any(
            pattern in label_lower
            for pattern in self.CYCLASE_LABEL_PATTERNS
        )

        if has_cyclase_pattern:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular lyase: cyclization/ring formation"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No intramolecular lyase (cyclase) pattern"
        )
