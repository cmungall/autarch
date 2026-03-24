"""Carbon-carbon bond hydrolase reaction classification.

C-C hydrolases (EC 3.7.1) cleave carbon-carbon bonds using water.
They typically act on diketones or keto-acids adjacent to carbonyl groups.

Pattern: R-CO-CH2-CO-R' + H2O → R-CO-OH + CH3-CO-R' (simplified)
"""

from abc import abstractmethod

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.molecules import (
    CHEBI_H2O,
)

# Common CC-hydrolase substrates
CHEBI_FUMARYLACETOACETATE = "CHEBI:18034"
CHEBI_ACETYLPYRUVATE = "CHEBI:15346"


class HydrolaseActingOnAcidCarbonCarbonBondsInKetonicSubstances(Hydrolase):
    """Legacy implementation for hydrolase acting on carbon-carbon bonds in ketonic substances.

    EC 3.7.1.x includes:
    - Fumarylacetoacetase: fumarylacetoacetate + H2O → acetoacetate + fumarate
    - Acetylpyruvate hydrolase: acetylpyruvate + H2O → acetate + pyruvate
    - Kynureninase: L-kynurenine + H2O → anthranilate + L-alanine

    These enzymes cleave C-C bonds adjacent to carbonyl groups.

    Examples:
    - Fumarylacetoacetate + H2O → acetoacetate + fumarate
    - Acetylpyruvate + H2O → acetate + pyruvate
    """

    EC_NUMBER_PREFIX = "3.7.1.-"  # C-C hydrolases

    @abstractmethod
    def _implementation_only(self) -> None:
        """Mark this legacy implementation class as abstract."""

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a C-C hydrolase.

        Strategy:
        1. Must use water as substrate
        2. Products should be two smaller carbon-containing fragments
        3. Look for characteristic diketone/keto-acid substrates
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must use water
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water substrate - C-C hydrolases are hydrolases"
            )

        # Check for characteristic substrates - ChEBI ID primary
        has_cc_substrate_chebi = any(
            p.chebi_id in {CHEBI_FUMARYLACETOACETATE, CHEBI_ACETYLPYRUVATE}
            for p in reaction.left_participants
        )

        # Fallback: label pattern matching
        cc_hydrolase_substrates = [
            "fumarylacetoacetate", "acetylpyruvate", "kynurenine",
            "oxaloacetate", "oxalosuccinate", "isocitrate",
            "dione", "diketone", "oxo"
        ]
        has_cc_substrate_label = any(s in substrate_str for s in cc_hydrolase_substrates)
        has_cc_substrate = has_cc_substrate_chebi or has_cc_substrate_label

        # Check for characteristic products (two carboxylic acids or ketones)
        carboxylic_products = [
            "acetate", "fumarate", "pyruvate", "oxalate",
            "succinate", "malate", "formate", "acetoacetate",
            "anthranilate", "alanine"
        ]
        product_matches = sum(1 for cp in carboxylic_products if cp in product_str)

        # Exclude reactions that release NH4+ (these are usually lyases or deaminases)
        has_nh4 = any(x in product_str for x in ["ammonium", "nh4", "ammonia"])

        if has_nh4:
            return ClassificationResult(
                is_member=False,
                explanation="Releases NH4+ - lyase or deaminase, not C-C hydrolase"
            )

        # Exclude thioesterases (release CoA)
        if "coa" in product_str:
            return ClassificationResult(
                is_member=False,
                explanation="Thioesterase (releases CoA) - not C-C hydrolase"
            )

        # Exclude esterases (produce alcohol + carboxylic acid from ester)
        ester_patterns = ["methyl ester", "ethyl ester", "acetoxy"]
        if any(ep in substrate_str for ep in ester_patterns):
            return ClassificationResult(
                is_member=False,
                explanation="Esterase - not C-C hydrolase"
            )

        # Exclude amidases (C-N cleavage, not C-C)
        amide_patterns = ["acetylaryl", "acetylamin", "amide"]
        if any(ap in substrate_str for ap in amide_patterns):
            return ClassificationResult(
                is_member=False,
                explanation="Amidase (C-N cleavage) - not C-C hydrolase"
            )

        # Exclude sulfur-containing reactions (H2S release)
        if "sulfide" in product_str or "h2s" in product_str:
            return ClassificationResult(
                is_member=False,
                explanation="Sulfur lyase - not C-C hydrolase"
            )

        # C-C hydrolysis typically produces two distinct fragments
        if has_cc_substrate and product_matches >= 2:
            # Find which substrate pattern matched for explanation
            matched_substrate = next(
                (s for s in cc_hydrolase_substrates if s in substrate_str),
                "diketone"
            )
            return ClassificationResult(
                is_member=True,
                explanation=f"C-C hydrolase: {matched_substrate} cleavage"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No C-C bond hydrolysis pattern"
        )
