"""Decarboxylase reaction classification using pattern DSL.

Decarboxylases are a specific type of lyase that remove CO2 from substrates.

History
-------

## 2025-12-21

Refactored to use SMARTS-based carboxylic acid detection (Moiety.CARBOXYL) via
Participant.is_carboxylic_acid() for substrate recognition. Name-based patterns
retained as supplementary detection.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_O2,
    co2,
    h_plus,
    p,
)


class CarboxyLyase(Lyase):
    """carboxy-lyase

    Examples:
    - Pyruvate decarboxylase: pyruvate → acetaldehyde + CO2
    - Amino acid decarboxylases: histidine → histamine + CO2
    - Ribulose bisphosphate carboxylase (RuBisCO): can work in reverse
    """

    GO_ID = "GO:0016831"  # carboxy-lyase activity
    EC_NUMBER_PREFIX = "4.1.1.-"

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple decarboxylation: substrate → product + CO2
        var("substrate") >> var("product") + p(co2) + optional(h_plus),
        # With multiple products: substrate → product1 + product2 + CO2
        var("substrate")
        >> var("product1") + var("product2") + p(co2) + optional(h_plus),
        # With cofactor: substrate + cofactor → product + CO2 + modified_cofactor
        var("substrate") + var("cofactor")
        >> var("product") + p(co2) + var("modified_cofactor") + optional(h_plus),
        # Multiple substrates: substrate1 + substrate2 → product + CO2
        var("substrate1") + var("substrate2")
        >> var("product") + p(co2) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a decarboxylase using pattern matching.

        Strategy:
        1. Must be a lyase (parent class)
        2. Must produce CO2 (not consume it)
        3. Must involve fragmentation
        4. Try pattern matching for specific decarboxylation patterns
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""

        # First check if it's a lyase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False, explanation=f"Not a lyase: {parent_result.explanation}"
            )

        # Must produce CO2 (on RIGHT side only)
        co2_chebi = CHEBI_CO2
        has_co2_product = any(
            participant.chebi_id == co2_chebi
            for participant in reaction.right_participants
        )

        if not has_co2_product:
            return ClassificationResult(
                is_member=False, explanation="No CO2 produced - not decarboxylation"
            )

        # Must NOT consume CO2 (it should be produced, not consumed)
        has_co2_reactant = any(
            participant.chebi_id == co2_chebi
            for participant in reaction.left_participants
        )

        if has_co2_reactant:
            return ClassificationResult(
                is_member=False,
                explanation="CO2 consumed rather than produced - likely carboxylation not decarboxylation",
            )

        # Exclude complex reactions that aren't simple decarboxylations
        
        # Exclude monooxygenase reactions (O2 + NAD(P)H involvement)
        has_oxygen = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )
        has_nadh_or_nadph = any(
            p.chebi_id in [CHEBI_NADPH, CHEBI_NADH]for p in reaction.left_participants
        )
        
        if has_oxygen and has_nadh_or_nadph:
            return ClassificationResult(
                is_member=False,
                explanation="Monooxygenase reaction - not simple decarboxylation"
            )

        # Exclude oxidative decarboxylation with O2 alone
        # Pattern: substrate + O2 → products + CO2 (tryptophan oxygenase, etc.)
        if has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="Oxidative decarboxylation using O2 - oxygenase not decarboxylase"
            )

        # Exclude oxidoreductases (cytochrome-coupled reactions)
        # Pattern: substrate + cytochrome → products + CO2 + reduced cytochrome
        cytochrome_patterns = ["cytochrome", "ferricytochrome", "ferrocytochrome"]
        has_cytochrome = any(pat in label_lower for pat in cytochrome_patterns)
        if has_cytochrome:
            return ClassificationResult(
                is_member=False,
                explanation="Cytochrome-coupled oxidoreductase - not decarboxylase"
            )

        # Exclude formate oxidation (formate dehydrogenase)
        CHEBI_FORMATE = "CHEBI:15740"
        has_formate = any(
            p.chebi_id == CHEBI_FORMATE
            for p in reaction.left_participants
        ) or "formate" in label_lower
        if has_formate:
            return ClassificationResult(
                is_member=False,
                explanation="Formate oxidation - oxidoreductase not decarboxylase"
            )

        # Exclude fatty acid/polyketide synthase (malonyl-ACP condensation)
        acp_patterns = ["[acp]", "malonyl-acp", "acyl-[acp]", "acyl carrier protein"]
        has_acp = any(pat in label_lower for pat in acp_patterns)
        if has_acp:
            return ClassificationResult(
                is_member=False,
                explanation="ACP-dependent condensation - fatty acid synthase not decarboxylase"
            )

        # Exclude dehydratase reactions that also release CO2
        # True decarboxylases: substrate → product + CO2
        # Dehydratases: substrate → product + CO2 + H2O (water elimination + CO2)
        has_water_product = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )
        if has_water_product:
            return ClassificationResult(
                is_member=False,
                explanation="Dehydratase with CO2 release - not simple decarboxylase"
            )

        diff = ReactionDiff(reaction)

        # Must involve fragmentation (molecule breaks apart releasing CO2)
        if not diff.is_fragmentation:
            # Sometimes molecule counts balance if a cofactor is involved
            if diff.n_product_molecules == diff.n_reactant_molecules:
                # Check if pattern matches with cofactor
                match = match_patterns(reaction, self.PATTERNS[2:3], strict=False)
                if not match.matched:
                    return ClassificationResult(
                        is_member=False,
                        explanation="No fragmentation - not typical decarboxylation",
                    )
            else:
                return ClassificationResult(
                    is_member=False,
                    explanation="No fragmentation - not typical decarboxylation",
                )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Calculate how many fragments are produced
        fragments_produced = diff.n_product_molecules - diff.n_reactant_molecules

        # Build explanation
        extra_info = []

        if fragments_produced == 1:
            extra_info.append("classic 1→2 decarboxylation")
        elif fragments_produced > 1:
            extra_info.append(
                f"multiple fragments ({diff.n_reactant_molecules}→{diff.n_product_molecules})"
            )

        # Check for specific substrates
        known_decarboxylation_substrates = {
            "CHEBI:15361": "pyruvate",
            "CHEBI:27432": "histidine",
            "CHEBI:16810": "2-oxoglutarate",
            "CHEBI:30031": "succinate",
        }

        for participant in reaction.left_participants:
            if participant.chebi_id in known_decarboxylation_substrates:
                extra_info.append(
                    f"{known_decarboxylation_substrates[participant.chebi_id]} substrate"
                )
                break
            # SMARTS-based carboxylic acid detection (structural)
            elif participant.is_carboxylic_acid():
                extra_info.append("carboxyl-containing substrate")
                break

        # Label-based fallback for carboxyl-containing substrates
        if not extra_info and any(
            term in label_lower
            for term in ["carboxyl", "carboxylic", "amino acid", "keto acid"]
        ):
            extra_info.append("carboxyl-containing substrate")

        explanation = "Decarboxylase: CO2 elimination"
        if extra_info:
            explanation += f", {'/'.join(extra_info)}"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
