"""Carbon-sulfur lyase reaction classification - refactored for inheritance.

EC 4.4: Carbon-sulfur lyases catalyze the cleavage of C-S bonds.
These enzymes break carbon-sulfur bonds in various substrates.

This classifier inherits from Lyase which handles:
- Exclusion of hydrolase patterns (water as reactant)
- Exclusion of oxidoreductase patterns (redox cofactors)
- Exclusion of transferase patterns (SAM, CoA, UDP/GDP)
- General lyase stoichiometry checks
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase

# Sulfur-containing product ChEBI IDs
CHEBI_H2S = "CHEBI:29919"  # hydrogen sulfide


class CarbonSulfurLyase(Lyase):
    """carbon sulfur lyase

    Inherits general lyase exclusions from parent.

    Specific C-S lyase patterns:
    - Cystathionine cleavage
    - Cysteine/homocysteine desulfhydrase (→ H2S)
    - SAM → MTA (not SAM → SAH which is methyltransferase)
    - Sulfolactate/cysteate → sulfite
    - Glutathione conjugate cleavage
    """

    GO_ID = "GO:0016846"  # carbon-sulfur lyase activity
    EC_NUMBER_PREFIX = "4.4.-.-"

    # C-S specific exclusions (beyond parent Lyase exclusions)
    CS_EXCLUSION_INDICATORS = [
        # Disulfide exchange (thioredoxin-like, not C-S cleavage)
        "disulfide", "dithiol",
        # Prenyl transferases (cysteine prenylation)
        "geranylgeranyl", "farnesyl",
        # Lipoyl modifications
        "lipoyl",
        # Enzyme cysteine modifications
        "cysteine-containing protein",
    ]

    # Positive pattern substrates
    CS_LYASE_SUBSTRATES = [
        "cystathionine",
        "lactoylglutathione",
        "hydroxymethylglutathione",
        "propiothetin",
        "coenzyme m",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a carbon-sulfur lyase.

        Strategy:
        1. Must pass parent Lyase checks (excludes hydrolases, oxidoreductases)
        2. Apply C-S specific exclusions
        3. Match positive C-S lyase patterns
        """
        # First check parent Lyase class
        # This handles: water reactant, redox cofactors, transferase cofactors
        parent_result = super().check_membership_impl(reaction)

        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""
        all_names_str = label_lower  # For compatibility with existing patterns
        left_names_str = substrate_str
        right_names_str = product_str

        # Apply C-S specific exclusions
        has_cs_exclusion = any(
            indicator in all_names_str
            for indicator in self.CS_EXCLUSION_INDICATORS
        )

        if has_cs_exclusion:
            return ClassificationResult(
                is_member=False,
                explanation="Excluded pattern (disulfide/prenyl/lipoyl) - not C-S lyase"
            )

        # Exclude SAM → SAH (methyltransferase, not C-S lyase)
        # C-S lyase: SAM → MTA (thioadenosine)
        has_sah = "adenosylhomocysteine" in right_names_str or "adenosyl-l-homocysteine" in right_names_str
        if has_sah:
            return ClassificationResult(
                is_member=False,
                explanation="Methyltransferase (SAM→SAH) - not C-S lyase"
            )

        # Exclude CoA transfer (CoA on both sides)
        has_coa_both = "coa" in left_names_str and "coa" in right_names_str
        if has_coa_both:
            return ClassificationResult(
                is_member=False,
                explanation="CoA transfer - not C-S lyase"
            )

        # Exclude ACP reactions
        if "acp]" in all_names_str or "-acp" in all_names_str:
            return ClassificationResult(
                is_member=False,
                explanation="ACP reaction - transferase, not C-S lyase"
            )

        # POSITIVE PATTERNS

        # 1. Direct substrate match
        has_cs_substrate = any(
            substrate in left_names_str
            for substrate in self.CS_LYASE_SUBSTRATES
        )

        if has_cs_substrate:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: C-S substrate detected"
            )

        # 2. Cysteine/homocysteine + H2S product (desulfhydrase)
        has_cysteine = "cysteine" in left_names_str or "homocysteine" in left_names_str
        has_h2s = any(
            p.chebi_id == CHEBI_H2S
            for p in reaction.right_participants
        ) or "hydrogen sulfide" in right_names_str

        if has_cysteine and has_h2s:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: desulfhydrase (→ H2S)"
            )

        # 3. SAM → MTA (thioadenosine) - true C-S lyase
        has_sam = "adenosyl-l-methionine" in left_names_str
        has_mta = "thioadenosine" in right_names_str or "methylthioadenosine" in right_names_str

        if has_sam and has_mta:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: SAM → MTA"
            )

        # 4. Sulfolactate/cysteate → sulfite
        has_sulfo_substrate = "sulfolactate" in left_names_str or "cysteate" in left_names_str
        has_sulfite = "sulfite" in right_names_str

        if has_sulfo_substrate and has_sulfite:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: → sulfite"
            )

        # 5. Leukotriene C → glutathione + leukotriene A
        has_ltc = "leukotriene c" in left_names_str
        has_gsh_product = "glutathione" in right_names_str

        if has_ltc and has_gsh_product:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: leukotriene cleavage"
            )

        # 6. Methionine → methanethiol
        has_methionine = "methionine" in left_names_str and "adenosyl" not in left_names_str
        has_methanethiol = "methanethiol" in right_names_str

        if has_methionine and has_methanethiol:
            return ClassificationResult(
                is_member=True,
                explanation="Carbon-sulfur lyase: methionine γ-lyase"
            )

        # If parent said it's a lyase but no C-S specific pattern, not C-S lyase
        if parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation="Lyase but no C-S bond cleavage pattern"
            )

        return ClassificationResult(
            is_member=False,
            explanation=f"Not a lyase: {parent_result.explanation}"
        )
