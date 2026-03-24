"""Dehydrogenase reaction classification - refactored to use hierarchical inheritance.

Dehydrogenases are EC 1.1.1: Oxidoreductases acting on CH-OH group with NAD(P)+ as acceptor.
They catalyze oxidation of alcohols to aldehydes/ketones using NAD+ or NADP+.

This classifier inherits from OxidoreductaseActingOnCHOHGroupOfDonors which handles:
- Exclusion of CH-CH patterns (→ EnoylReductase EC 1.3)
- Exclusion of CH-NH2 patterns (→ AminoAcidDehydrogenase EC 1.4)
- Exclusion of electron carrier reactions
- Exclusion of disulfide reactions
- Exclusion of generic templates
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase_acting_on_ch_oh_group_of_donors import (
    OxidoreductaseActingOnCHOHGroupOfDonors,
)
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    h_plus,
    lactate,
    nad_plus,
    nadh,
    nadp_plus,
    nadph,
    p,
    pyruvate,
)


class OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor(
    OxidoreductaseActingOnCHOHGroupOfDonors
):
    """oxidoreductase acting on the CH-OH group of donors, NAD or NADP as acceptor

    EC 1.1.1.x includes:
    - Alcohol dehydrogenase: ethanol + NAD+ → acetaldehyde + NADH + H+
    - Lactate dehydrogenase: lactate + NAD+ → pyruvate + NADH + H+
    - Glucose-6-phosphate dehydrogenase: G6P + NADP+ → 6PG + NADPH + H+

    Key distinctions (now handled by parent OxidoreductaseActingOnCHOHGroupOfDonors):
    - NOT EC 1.1.3 (uses O2) → AlcoholOxidase
    - NOT EC 1.3 (CH-CH group) → EnoylReductase
    - NOT EC 1.4 (CH-NH2 group) → AminoAcidDehydrogenase
    """

    GO_ID = "GO:0016616"  # oxidoreductase activity, acting on CH-OH with NAD/NADP
    EC_NUMBER_PREFIX = "1.1.1.-"

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # NAD+ dependent oxidation: substrate + NAD+ → product + NADH + H+
        var("substrate") + p(nad_plus) >> var("product") + p(nadh) + p(h_plus),
        # NADP+ dependent oxidation: substrate + NADP+ → product + NADPH + H+
        var("substrate") + p(nadp_plus) >> var("product") + p(nadph) + p(h_plus),
        # Reverse (reduction) with NADH: substrate + NADH + H+ → product + NAD+
        var("substrate") + p(nadh) + p(h_plus) >> var("product") + p(nad_plus),
        # Reverse (reduction) with NADPH: substrate + NADPH + H+ → product + NADP+
        var("substrate") + p(nadph) + p(h_plus) >> var("product") + p(nadp_plus),
        # Without explicit H+: substrate + NAD+ → product + NADH
        var("substrate") + p(nad_plus) >> var("product") + p(nadh),
        # Lactate dehydrogenase specific
        p(lactate) + p(nad_plus) >> p(pyruvate) + p(nadh) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a dehydrogenase (EC 1.1.1).

        Strategy:
        1. Must pass OxidoreductaseActingOnCHOHGroupOfDonors checks (excludes CH-CH, CH-NH2, etc.)
        2. Must have proper NAD+→NADH or NADP+→NADPH conversion
        3. Must NOT use FAD (that's EC 1.1.99)
        4. Must NOT use O2 as acceptor (that's EC 1.1.3)
        """
        # Check parent class (handles CH-OH vs other substrate types)
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        # Must NOT use FAD (EC 1.1.99 uses other acceptors)
        fad_system = {CHEBI_FAD, CHEBI_FADH2}
        has_fad = any(
            p.chebi_id in fad_system
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_fad:
            return ClassificationResult(
                is_member=False,
                explanation="FAD-dependent - EC 1.1.99, not EC 1.1.1"
            )

        # Must NOT use O2 as acceptor (EC 1.1.3 = alcohol oxidases)
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )

        if has_o2:
            return ClassificationResult(
                is_member=False,
                explanation="O2 as acceptor - alcohol oxidase (EC 1.1.3), not dehydrogenase"
            )

        # Check the direction of NAD/NADP conversion
        # Forward: NAD(P)+ on left, NAD(P)H on right
        has_oxidized_left = any(
            p.chebi_id in [CHEBI_NAD_PLUS, CHEBI_NADP_PLUS]
            for p in reaction.left_participants
        )
        has_reduced_right = any(
            p.chebi_id in [CHEBI_NADH, CHEBI_NADPH]
            for p in reaction.right_participants
        )

        # Reverse: NAD(P)H on left, NAD(P)+ on right
        has_reduced_left = any(
            p.chebi_id in [CHEBI_NADH, CHEBI_NADPH]
            for p in reaction.left_participants
        )
        has_oxidized_right = any(
            p.chebi_id in [CHEBI_NAD_PLUS, CHEBI_NADP_PLUS]
            for p in reaction.right_participants
        )

        # Try pattern matching
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        if has_oxidized_left and has_reduced_right:
            # Forward direction: oxidation
            cofactor = (
                "NAD+"
                if any(p.chebi_id == CHEBI_NAD_PLUS for p in reaction.left_participants)
                else "NADP+"
            )
            explanation = f"Dehydrogenase: {cofactor}-dependent oxidation"
            if match.matched:
                explanation += " [pattern matched]"

            # Add substrate-specific notes
            explanation = self._add_substrate_notes(reaction, explanation)

            return ClassificationResult(is_member=True, explanation=explanation)

        if has_reduced_left and has_oxidized_right:
            # Reverse direction: reduction
            cofactor = (
                "NADH"
                if any(p.chebi_id == CHEBI_NADH for p in reaction.left_participants)
                else "NADPH"
            )
            explanation = f"Dehydrogenase: {cofactor}-dependent reduction (reverse)"
            if match.matched:
                explanation += " [pattern matched]"

            return ClassificationResult(is_member=True, explanation=explanation)

        return ClassificationResult(
            is_member=False,
            explanation="No clear NAD(P)+/NAD(P)H conversion pattern"
        )

    def _add_substrate_notes(self, reaction: Reaction, explanation: str) -> str:
        """Add substrate-specific notes to explanation."""
        # Known substrates
        substrate_notes = {
            "CHEBI:24996": ", lactate dehydrogenase",
            "CHEBI:16236": ", alcohol dehydrogenase",
            "CHEBI:17665": ", G6P dehydrogenase",
            "CHEBI:16908": ", isocitrate dehydrogenase",
            "CHEBI:16761": ", malate dehydrogenase",
        }

        for chebi_id, note in substrate_notes.items():
            if any(p.chebi_id == chebi_id for p in reaction.left_participants):
                explanation += note
                break

        return explanation
