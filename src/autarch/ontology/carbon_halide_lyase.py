"""Carbon-halide lyase reaction classification.

EC 4.5: Carbon-halide lyases catalyze the cleavage of C-halide bonds.
These enzymes remove halide ions (Cl-, Br-, F-, I-) from substrates
without hydrolysis.

Pattern: R-X → R' + HX (where X = halide)
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase


class CarbonHalideLyase(Lyase):
    """carbon-halide lyase

    EC 4.5.x.x includes:
    - DDT-dehydrochlorinase (4.5.1.1)
    - 3-chloro-D-alanine dehydrochlorinase (4.5.1.2)
    - Dichloromethane dehalogenase (4.5.1.3)
    - Haloalkane dehalogenase-like lyases

    Examples:
    - DDT = DDE + HCl
    - 3-chloro-D-alanine = pyruvate + NH3 + HCl
    """

    GO_ID = "GO:0016848"  # carbon-halide lyase activity
    EC_NUMBER_PREFIX = "4.5.-.-"

    # Halide-related ChEBI IDs
    CHEBI_CHLORIDE = "CHEBI:17996"
    CHEBI_BROMIDE = "CHEBI:15858"
    CHEBI_FLUORIDE = "CHEBI:17051"
    CHEBI_IODIDE = "CHEBI:16382"
    CHEBI_HCL = "CHEBI:17883"

    HALIDE_CHEBI_IDS = {
        CHEBI_CHLORIDE, CHEBI_BROMIDE, CHEBI_FLUORIDE, CHEBI_IODIDE, CHEBI_HCL,
    }

    HALIDE_LABEL_PATTERNS = [
        "dehalogenase", "dehydrochlorinase", "dechlorinase",
        "haloalkane", "haloacid",
        "chloro", "bromo", "fluoro", "iodo",
        "ddt", "dichloromethane", "halide",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a carbon-halide lyase.

        Strategy:
        1. Must pass parent Lyase checks
        2. Look for halide release (products)
        3. Look for label patterns indicating dehalogenation
        """
        parent_result = super().check_membership_impl(reaction)

        label_lower = reaction.label.lower() if reaction.label else ""

        # Check for halide as product (ChEBI-based)
        has_halide_product = any(
            p.chebi_id in self.HALIDE_CHEBI_IDS
            for p in reaction.right_participants
        )

        # Check for label-based halide patterns
        has_halide_label = any(
            pattern in label_lower
            for pattern in self.HALIDE_LABEL_PATTERNS
        )

        if has_halide_product or has_halide_label:
            if parent_result.is_member or has_halide_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="Carbon-halide lyase: halide release detected"
                )

        if parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation="Lyase but no C-halide cleavage pattern"
            )

        return ClassificationResult(
            is_member=False,
            explanation=f"Not a lyase: {parent_result.explanation}"
        )
