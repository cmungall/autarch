"""Monooxygenase reaction classification using pattern DSL."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    h_plus,
    p,
    water,
)

# Iron redox for oxidase exclusion
CHEBI_FE2 = "CHEBI:29033"  # Fe(II)
CHEBI_FE3 = "CHEBI:29034"  # Fe(III)

# Ferredoxin for desaturase exclusion
CHEBI_FERREDOXIN_OX = "CHEBI:33737"
CHEBI_FERREDOXIN_RED = "CHEBI:33738"

# Quinone/hydroquinone for oxidase exclusion
CHEBI_HYDROQUINONE = "CHEBI:17594"
CHEBI_SEMIQUINONE = "CHEBI:17977"


class Monooxygenase(Oxidoreductase):
    """monooxygenase

    Examples:
    - Cytochrome P450: substrate + O2 + NADPH + H+ → hydroxylated substrate + H2O + NADP+
    - Catechol dioxygenase: catechol + O2 → muconic acid
    - Tyrosinase: tyrosine + O2 → DOPA + H2O
    - Lipoxygenase: arachidonic acid + O2 → hydroperoxyeicosatetraenoic acid
    """

    GO_ID = "GO:0004497"  # monooxygenase activity
    EC_NUMBER_PREFIX = None  # Broad GO xref only; do not inherit Oxidoreductase exact EC prefix
    EC_BROAD_XREFS = ["1.-.-.-"]  # GO xref is broadMatch, not exact EC equivalence

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Monooxygenase: substrate + O2 + electron donor → hydroxylated + H2O + oxidized donor
        var("substrate") + var("oxygen") + var("electron_donor") + optional(h_plus)
        >> var("hydroxylated_product") + p(water) + var("oxidized_donor"),
        # Dioxygenase: substrate + O2 → oxygenated product (both O atoms incorporated)
        var("substrate") + var("oxygen") >> var("oxygenated_product"),
        # P450-type: substrate + O2 + NADPH + H+ → hydroxylated + H2O + NADP+
        var("substrate") + var("o2") + var("nadph") + p(h_plus)
        >> var("hydroxylated") + p(water) + var("nadp_plus"),
        # Simple oxygenation: organic substrate + O2 → oxidized substrate + H2O
        var("organic_substrate") + var("molecular_oxygen")
        >> var("oxidized_substrate") + p(water) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxygenase using pattern matching.

        Strategy:
        1. Must be an oxidoreductase (parent class)
        2. Must use molecular oxygen (O2) as oxidant
        3. Must incorporate oxygen into substrate (not just remove H)
        4. Look for characteristic oxygenation patterns
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Must use molecular oxygen as oxidant
        oxygen_indicators = {
            CHEBI_O2,  # dioxygen (O2)
        }

        has_oxygen = any(
            p.chebi_id in oxygen_indicators
            for p in reaction.left_participants
        )

        if not has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="No molecular oxygen - not monooxygenase"
            )

        # Look for oxygen incorporation products (hydroxylation, oxygenation)

        has_oxygenated_product = False

        # Alternative: check for water production (common in monooxygenases)
        has_water_product = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )

        if not (has_oxygenated_product or has_water_product):
            return ClassificationResult(
                is_member=False,
                explanation="No oxygen incorporation pattern detected"
            )

        # Apply exclusions before classifying as oxygenase

        # Exclude simple oxidases (don't incorporate oxygen into substrate)
        # These typically produce H2O2 instead of incorporating O
        has_h2o2_product = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        )
        
        if has_h2o2_product and not has_oxygenated_product:
            return ClassificationResult(
                is_member=False,
                explanation="Oxidase producing H2O2 - not monooxygenase (no O incorporation)"
            )

        # Exclude simple dehydrogenases (shouldn't use O2 as primary oxidant)
        nad_system = {
            CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH  # NAD+, NADH, NADP+, NADPH
        }
        
        has_nad_primary = any(
            p.chebi_id in nad_system
            for p in reaction.left_participants
        ) and not any(
            p.chebi_id in nad_system
            for p in reaction.right_participants
        )
        
        # Allow NAD(P)H as electron donor for monooxygenases, but not as primary substrate
        if has_nad_primary:
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P)-dependent dehydrogenase - not oxygenase"
            )

        # Exclude peroxidases (use H2O2, not O2)
        has_h2o2_substrate = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.left_participants
        )

        if has_h2o2_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Peroxidase (uses H2O2) - not monooxygenase (uses O2)"
            )

        # Exclude iron oxidases (Fe2+ + O2 → Fe3+ + H2O)
        has_fe2_reactant = any(
            p.chebi_id == CHEBI_FE2
            for p in reaction.left_participants
        )
        has_fe3_product = any(
            p.chebi_id == CHEBI_FE3
            for p in reaction.right_participants
        )
        if has_fe2_reactant and has_fe3_product:
            return ClassificationResult(
                is_member=False,
                explanation="Iron oxidase (Fe2+ → Fe3+) - O2 reduced to H2O, not incorporated"
            )

        # Exclude ferredoxin-dependent desaturases (O2 → H2O, no O incorporation)
        has_ferredoxin_red = any(
            p.chebi_id == CHEBI_FERREDOXIN_RED
            for p in reaction.left_participants
        )
        has_ferredoxin_ox = any(
            p.chebi_id == CHEBI_FERREDOXIN_OX
            for p in reaction.right_participants
        )
        if has_ferredoxin_red and has_ferredoxin_ox:
            return ClassificationResult(
                is_member=False,
                explanation="Ferredoxin-dependent desaturase - O2 reduced to H2O, not incorporated"
            )

        # Exclude quinone oxidases (hydroquinone + O2 → semiquinone + H2O)
        has_hydroquinone = any(
            p.chebi_id == CHEBI_HYDROQUINONE
            for p in reaction.left_participants
        )
        has_semiquinone = any(
            p.chebi_id == CHEBI_SEMIQUINONE
            for p in reaction.right_participants
        )
        if has_hydroquinone and has_semiquinone:
            return ClassificationResult(
                is_member=False,
                explanation="Quinone oxidase - O2 reduced to H2O, not incorporated"
            )

        # Exclude 2-oxoglutarate-dependent dioxygenases (EC 1.14.11.x)
        # These are dioxygenases (GO:0016706), not monooxygenases (GO:0004497)
        # Pattern: substrate + 2-oxoglutarate + O2 → product + succinate + CO2
        CHEBI_2_OXOGLUTARATE = "CHEBI:16810"  # 2-oxoglutarate/alpha-ketoglutarate
        CHEBI_SUCCINATE = "CHEBI:15741"  # succinate
        CHEBI_CO2 = "CHEBI:16526"  # carbon dioxide
        has_2_oxoglutarate = any(
            p.chebi_id == CHEBI_2_OXOGLUTARATE
            for p in reaction.left_participants
        ) or "oxoglutarate" in substrate_str
        has_succinate = any(
            p.chebi_id == CHEBI_SUCCINATE
            for p in reaction.right_participants
        ) or "succinate" in product_str
        has_co2 = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.right_participants
        ) or "co2" in product_str or "carbon dioxide" in product_str
        if has_2_oxoglutarate and has_succinate and has_co2:
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate-dependent dioxygenase - not monooxygenase"
            )

        # Exclude bilirubin oxidases (EC 1.3.3.5) - oxidizes bilirubin without O incorporation
        bilirubin_patterns = ["bilirubin", "biliverdin"]
        has_bilirubin = any(pat in label_lower for pat in bilirubin_patterns)
        if has_bilirubin:
            return ClassificationResult(
                is_member=False,
                explanation="Bilirubin oxidase - electron transfer not O incorporation"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Determine specific type
        explanation = "Monooxygenase: oxygen incorporation"
        
        # Check for electron donors to distinguish mono- vs dioxygenase
        has_electron_donor = any(
            p.chebi_id in {CHEBI_NADPH, CHEBI_NADH} for p in reaction.left_participants
        )
        
        if has_electron_donor and has_water_product:
            explanation += " (monooxygenase)"
        elif not has_electron_donor and not has_water_product:
            explanation += " (dioxygenase)"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
