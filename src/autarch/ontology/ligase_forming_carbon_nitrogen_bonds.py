"""Synthase reaction classification using pattern DSL.

Synthases are ligases that form bonds using ATP energy.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ligase import Ligase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_COA,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    adp,
    ammonia,
    ammonium,
    amp,
    atp,
    co2,
    diphosphate,
    h_plus,
    p,
    phosphate,
    water,
)

# tRNA polymer for aminoacyl-tRNA synthetase exclusion
CHEBI_TRNA = "CHEBI:78442"


class LigaseFormingCarbonNitrogenBonds(Ligase):
    """ligase forming carbon-nitrogen bonds

    Examples:
    - Glutamine synthetase: glutamate + NH3 + ATP → glutamine + ADP + Pi
    - Asparagine synthetase: aspartate + glutamine + ATP → asparagine + glutamate + AMP + PPi
    - Carbamoyl phosphate synthetase: NH3 + CO2 + 2ATP → carbamoyl phosphate + 2ADP + Pi

    Note: GO:0016879 specifically refers to C-N bond formation, not general synthesis
    """

    GO_ID = "GO:0016879"  # ligase activity, forming carbon-nitrogen bonds
    EC_NUMBER_PREFIX = "6.3.-.-"  # Ligases forming carbon-nitrogen bonds

    # Define patterns using DSL with operator overloading
    # Type as list[Reaction] for mypy compatibility
    PATTERNS: list[Reaction] = [
        # General ATP-dependent synthesis: A + B + ATP → A-B + ADP + Pi + H+?
        var("substrate1") + var("substrate2") + p(atp)
        >> var("product") + p(adp) + p(phosphate) + optional(h_plus),
        # With water involvement: A + B + ATP + H2O → A-B + ADP + Pi + H+?
        var("substrate1") + var("substrate2") + p(atp) + p(water)
        >> var("product") + p(adp) + p(phosphate) + optional(h_plus),
        # AMP + PPi variant: A + B + ATP → A-B + AMP + PPi + H+?
        var("substrate1") + var("substrate2") + p(atp)
        >> var("product") + p(amp) + p(diphosphate) + optional(h_plus),
        # CO2 fixation: substrate + CO2 + ATP → product + ADP + Pi + H+?
        var("substrate") + p(co2) + p(atp)
        >> var("product") + p(adp) + p(phosphate) + optional(h_plus),
        # Amino acid synthesis: substrate + NH3 + ATP → product + ADP + Pi + H+?
        var("substrate") + p(ammonia) + p(atp)
        >> var("product") + p(adp) + p(phosphate) + optional(h_plus),
        # Ammonium variant: substrate + NH4+ + ATP → product + ADP + Pi
        var("substrate") + p(ammonium) + p(atp)
        >> var("product") + p(adp) + p(phosphate),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a synthase using pattern matching.

        Strategy:
        1. First check if it's a ligase (ATP-dependent bond formation)
        2. Strict validation of ATP requirement and fusion pattern
        3. Try pattern matching for specific synthesis patterns
        4. Apply procedural logic for complex validation
        """
        # First check if it's a ligase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a ligase: {parent_result.explanation}",
            )

        # Must have ATP as reactant (LEFT side only)
        atp_in_reactants = any(
            p.chebi_id == CHEBI_ATP for p in reaction.left_participants
        )
        if not atp_in_reactants:
            return ClassificationResult(
                is_member=False,
                explanation="No ATP consumption - synthases are ATP-dependent ligases",
            )

        # Must NOT have ATP in products (it should be consumed)
        atp_in_products = any(
            p.chebi_id == CHEBI_ATP for p in reaction.right_participants
        )
        if atp_in_products:
            return ClassificationResult(
                is_member=False, explanation="ATP not consumed - appears in products"
            )

        # Must produce ADP or AMP (not just transfer ATP)
        has_adp_product = any(
            p.chebi_id == CHEBI_ADP for p in reaction.right_participants
        )
        has_amp_product = any(
            p.chebi_id == CHEBI_AMP for p in reaction.right_participants
        )

        if not (has_adp_product or has_amp_product):
            return ClassificationResult(
                is_member=False,
                explanation="No ADP/AMP produced - not ATP hydrolysis/utilization",
            )

        # === C-N SPECIFIC EXCLUSIONS ===

        # 1. Exclude aminoacyl-tRNA synthetases (form C-O ester bonds, not C-N)
        has_trna = any(
            p.chebi_id == CHEBI_TRNA
            for p in reaction.left_participants + reaction.right_participants
        )
        if has_trna:
            return ClassificationResult(
                is_member=False,
                explanation="tRNA synthetase - forms C-O ester bond, not C-N"
            )

        # 2. Exclude CoA ligases (form C-S thioester bonds, not C-N)
        has_coa_reactant = any(
            p.chebi_id == CHEBI_COA
            for p in reaction.left_participants
        )
        if has_coa_reactant:
            return ClassificationResult(
                is_member=False,
                explanation="CoA ligase - forms C-S thioester bond, not C-N"
            )

        # 3. Exclude oxygenases/oxidases (use O2)
        has_oxygen = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )
        if has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygenase/oxidase - uses O2, not C-N ligase"
            )

        # Check for nitrogen source (C-N ligases should involve nitrogen)
        has_nitrogen_source = any(
            s.chebi_id in [CHEBI_NH3, CHEBI_NH4]  # NH3, NH4+
            for s in reaction.left_participants
        )
        amino_acid_donor_terms = {
            "alanine",
            "arginine",
            "asparagine",
            "aspartate",
            "cysteine",
            "glutamate",
            "glutamine",
            "glycine",
            "histidine",
            "isoleucine",
            "leucine",
            "lysine",
            "methionine",
            "phenylalanine",
            "proline",
            "serine",
            "threonine",
            "tryptophan",
            "tyrosine",
            "valine",
        }
        has_amino_acid_donor = any(
            participant.name
            and any(term in participant.name.lower() for term in amino_acid_donor_terms)
            for participant in reaction.left_participants
        )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        if match.matched:
            # Pattern matched - determine synthesis type
            if "substrate1" in match.bindings and "substrate2" in match.bindings:
                s1 = match.bindings.get("substrate1", "substrate1")
                s2 = match.bindings.get("substrate2", "substrate2")
                explanation = f"Synthase: ATP-dependent synthesis of {s1} + {s2}"
            elif "substrate" in match.bindings:
                substrate = match.bindings.get("substrate", "substrate")

                # Check for specific synthesis types based on co-reactants
                if match.has_unmatched(CHEBI_CO2):  # CO2
                    explanation = (
                        f"Synthase: ATP-dependent CO2 fixation with {substrate}"
                    )
                elif match.has_unmatched(CHEBI_NH3) or match.has_unmatched(
                    CHEBI_NH4
                ):  # NH3/NH4+
                    explanation = (
                        f"Synthase: ATP-dependent amino acid synthesis from {substrate}"
                    )
                else:
                    explanation = f"Synthase: ATP-dependent synthesis with {substrate}"
            else:
                explanation = "Synthase: ATP-dependent bond formation"

            # Add product info if available
            if has_adp_product:
                explanation += " (ADP + Pi)"
            elif has_amp_product:
                explanation += " (AMP + PPi)"

            return ClassificationResult(is_member=True, explanation=explanation)

        # No pattern matched - fall back to procedural checks
        diff = ReactionDiff(reaction)

        # Should involve fusion (combining substrates)
        if not diff.is_fusion:
            return ClassificationResult(
                is_member=False, explanation="No fusion - synthases combine substrates"
            )

        # Check for common synthesis patterns

        # Pattern 1: CO2 fixation (carbon-carbon bond formation)
        has_co2_reactant = any(
            s.chebi_id == CHEBI_CO2 for s in reaction.left_participants
        )

        if has_co2_reactant:
            return ClassificationResult(
                is_member=True, explanation="Synthase: ATP-dependent CO2 fixation"
            )

        # Pattern 2: Amino acid synthesis (C-N bond formation)  
        # has_nitrogen_source already defined above
        if (has_nitrogen_source or has_amino_acid_donor) and "N" in diff.product_elements:
            return ClassificationResult(
                is_member=True,
                explanation="Synthase: ATP-dependent amino acid/amide synthesis",
            )

        # Pattern 3: General ATP-dependent synthesis
        # Multiple substrates + ATP → single product + ADP + Pi
        if (
            diff.n_reactant_molecules >= 3
            and diff.n_product_molecules <= diff.n_reactant_molecules - 1
        ):
            return ClassificationResult(
                is_member=True,
                explanation=f"Synthase: ATP-dependent synthesis ({diff.n_reactant_molecules}→{diff.n_product_molecules})",
            )

        # Must involve nitrogen incorporation for C-N bond formation
        if not (has_nitrogen_source or has_amino_acid_donor) and "N" not in diff.product_elements:
            return ClassificationResult(
                is_member=False,
                explanation="No nitrogen source - C-N ligases require nitrogen incorporation"
            )
        
        # If it passes all checks, it's a C-N ligase
        return ClassificationResult(
            is_member=True,
            explanation="Carbon-Nitrogen Ligase: ATP-dependent C-N bond formation with nitrogen incorporation",
        )
