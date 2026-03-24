"""Translocase reaction classification using pattern DSL.

Translocases catalyze the movement of ions or molecules across membranes.
EC 7.x.x.x in the enzyme classification system (added in 2018).

NOTE: Current evaluation dataset contains ZERO translocase reactions (GO:0015085).
While ATP-driven transporters exist in the data, they're not classified as 
"translocase activity" - likely a GO term annotation issue.

Previous implementation had 8 false positives with 0 true positives.
Disabled until better test data available.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    adp,
    amp,
    atp,
    diphosphate,
    h_plus,
    p,
    phosphate,
    water,
)


class Translocase(ReactionClass):
    """Translocase reaction classifier using pattern DSL.

    Translocases catalyze transmembrane movement:
    A(out) + ATP → A(in) + ADP + Pi

    Common examples:
    - Ion pumps: Na+/K+-ATPase, Ca2+-ATPase
    - ABC transporters: drug efflux pumps
    - Proton pumps: H+-ATPase
    - Symporters/antiporters (if coupled to ATP)

    Note: These often look like ATP hydrolysis but the key is
    the (in)/(out) notation indicating membrane transport.
    """

    GO_ID = "GO:0015085"  # transporter activity
    EC_NUMBER_PREFIX = "7.-.-.-"

    # Define patterns using DSL with operator overloading
    # Type as list[Reaction] for mypy compatibility
    PATTERNS: list[Reaction] = [
        # Basic ATP-driven transport: substrate + ATP + H2O → product + ADP + Pi + H+?
        var("substrate") + p(atp) + p(water)
        >> var("product") + p(adp) + p(phosphate) + optional(h_plus),
        # Ion pump patterns (same ion on both sides with location change)
        # Na+ transport: Na+(out) + ATP + H2O → Na+(in) + ADP + Pi
        var("ion_out") + p(atp) + p(water) >> var("ion_in") + p(adp) + p(phosphate),
        # Multi-ion transport with protons: ion + ATP + H2O → ion + ADP + Pi + H+
        var("ion1") + p(atp) + p(water)
        >> var("ion2") + p(adp) + p(phosphate) + p(h_plus),
        # AMP + PPi variant: substrate + ATP → product + AMP + PPi
        var("substrate") + p(atp) >> var("product") + p(amp) + p(diphosphate),
        # Simple transport without explicit water: substrate + ATP → product + ADP + Pi
        var("substrate") + p(atp) >> var("product") + p(adp) + p(phosphate),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a translocase using pattern matching.

        Strategy:
        1. First check if it's a transport reaction using location information
        2. Strict validation of ATP requirement for active transport
        3. Try pattern matching for specific transport patterns
        4. Apply procedural logic for complex validation
        """
        # First check if this is a transport reaction with location information
        if reaction.is_transport_reaction():
            # Check if ATP is involved (energy-dependent transport)
            has_atp = any(
                p.chebi_id == CHEBI_ATP for p in reaction.left_participants
            )
            has_adp = any(
                p.chebi_id == CHEBI_ADP for p in reaction.right_participants
            )

            if has_atp and has_adp:
                transported = reaction.get_transported_molecules()
                if transported:
                    mol_names = [t[0] for t in transported]
                    return ClassificationResult(
                        is_member=True,
                        explanation=f"Translocase: ATP-driven transport of {', '.join(mol_names[:2])}",
                    )

            # Some translocases don't use ATP (e.g., facilitated diffusion)
            # but we're focusing on active transport for now
            return ClassificationResult(
                is_member=False, explanation="Transport reaction but no ATP involvement"
            )

        # Must have ATP as reactant (most translocases are ATP-driven)
        atp_in_reactants = any(
            p.chebi_id == CHEBI_ATP for p in reaction.left_participants
        )
        if not atp_in_reactants:
            return ClassificationResult(
                is_member=False,
                explanation="No ATP - most translocases require ATP for active transport",
            )

        # Must NOT have ATP in products (it should be consumed)
        atp_in_products = any(
            p.chebi_id == CHEBI_ATP for p in reaction.right_participants
        )
        if atp_in_products:
            return ClassificationResult(
                is_member=False,
                explanation="ATP not consumed - not ATP-driven transport",
            )

        # Must produce ADP or AMP (energy coupling)
        has_adp_product = any(
            p.chebi_id == CHEBI_ADP for p in reaction.right_participants
        )
        has_amp_product = any(
            p.chebi_id == CHEBI_AMP for p in reaction.right_participants
        )

        if not (has_adp_product or has_amp_product):
            return ClassificationResult(
                is_member=False,
                explanation="No ADP/AMP produced - not ATP hydrolysis coupling",
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        if match.matched:
            # Check if substrate and product could represent the same transported molecule
            if "substrate" in match.bindings and "product" in match.bindings:
                substrate_id = match.bindings["substrate"]
                product_id = match.bindings["product"]

                # If same CHEBI ID, it's likely transport (same molecule, different location)
                if substrate_id == product_id:
                    return ClassificationResult(
                        is_member=True,
                        explanation=f"Translocase: ATP-driven transport of {substrate_id}",
                    )
                else:
                    # Different molecules - check if they're related ions or states
                    return ClassificationResult(
                        is_member=True,
                        explanation=f"Translocase: ATP-driven transport/transformation from {substrate_id} to {product_id}",
                    )

            # Check for ion transport patterns
            elif "ion_out" in match.bindings and "ion_in" in match.bindings:
                ion_out = match.bindings["ion_out"]
                ion_in = match.bindings["ion_in"]

                if ion_out == ion_in:
                    return ClassificationResult(
                        is_member=True,
                        explanation=f"Translocase: ATP-driven ion transport ({ion_out})",
                    )
                else:
                    return ClassificationResult(
                        is_member=True,
                        explanation=f"Translocase: ATP-driven transport ({ion_out} → {ion_in})",
                    )

            # General pattern match
            explanation = "Translocase: ATP-driven transport"
            if has_adp_product:
                explanation += " (ADP + Pi)"
            elif has_amp_product:
                explanation += " (AMP + PPi)"

            return ClassificationResult(is_member=True, explanation=explanation)

        # No pattern matched - fall back to procedural checks
        diff = ReactionDiff(reaction)

        # Get all CHEBI IDs on each side
        left_chebis = {s.chebi_id for s in reaction.left_participants if s.chebi_id}
        right_chebis = {s.chebi_id for s in reaction.right_participants if s.chebi_id}

        # Exclude ATP/ADP/AMP and common cofactors
        energy_molecules = {
            CHEBI_ATP,  # ATP
            CHEBI_ADP,  # ADP
            CHEBI_AMP,  # AMP
            CHEBI_PHOSPHATE,  # phosphate
            CHEBI_DIPHOSPHATE,  # diphosphate
            CHEBI_H2O,  # water
            CHEBI_H_PLUS,  # H+
        }

        # Remove energy molecules from comparison
        transported_left = left_chebis - energy_molecules
        transported_right = right_chebis - energy_molecules

        # Check if same molecule appears on both sides (suggests transport)
        common_molecules = transported_left & transported_right

        if common_molecules:
            # Same non-energy molecule on both sides with ATP consumption
            # This is likely transport across a membrane
            return ClassificationResult(
                is_member=True,
                explanation="Translocase: ATP-driven transport (same molecule both sides)",
            )

        # Pattern: ATP + H2O + X(out) → ADP + Pi + X(in) + H+
        # This is classic ATP-driven transport
        if diff.has_water_reactant and diff.is_fragmentation:
            # Check if it's ATP hydrolysis coupled to transport
            has_atp_to_adp = (
                any(s.chebi_id == CHEBI_ATP for s in reaction.left_participants)
                and any(
                    s.chebi_id == CHEBI_ADP for s in reaction.right_participants
                )
                and any(
                    s.chebi_id == CHEBI_PHOSPHATE for s in reaction.right_participants
                )
            )

            if has_atp_to_adp:
                # ATP hydrolysis is occurring
                # If there are other molecules involved, might be transport
                if len(transported_left) > 0 or len(transported_right) > 0:
                    return ClassificationResult(
                        is_member=True, explanation="Translocase: ATP-driven transport"
                    )

        return ClassificationResult(
            is_member=False, explanation="No clear transport pattern detected"
        )
