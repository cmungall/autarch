"""Simplified Hydrolase reaction classification.

Hydrolases break bonds using water: A-B + H2O → A-OH + B-H

History
-------

## 2025-12-22 (v2)

Fixed water detection to support multiple SMILES encodings:
- "O", "[H]O[H]", "[OH2]", "O([H])[H]"
Previously only matched CHEBI_H2O or SMILES "O", missing some valid reactions.

## 2025-12-22 (v1)

Fixed inheritance bug: check_membership_impl() now uses Hydrolase.PATTERNS
explicitly instead of self.PATTERNS. Previously, subclasses calling super()
would incorrectly use their own patterns for the parent check, causing valid
hydrolysis reactions to be rejected.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class Hydrolase(ReactionClass):
    """hydrolase"""
    
    GO_ID = "GO:0016787"  # hydrolase activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "3.-.-.-"  # All hydrolases
    
    # Ultra-simple patterns
    PATTERNS: list[Reaction] = [
        # Basic hydrolysis: substrate + H2O → product1 + product2 + H+?
        var("substrate") + p(water) >> var("product1") + var("product2") + optional(h_plus),
        # Without H+
        var("substrate") + p(water) >> var("product1") + var("product2"),
        # Three products
        var("substrate") + p(water) >> var("product1") + var("product2") + var("product3"),
    ]
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a hydrolase using ONLY strict pattern matching.
        
        SIMPLE, DECLARATIVE: Only reactions that match hydrolysis patterns.
        No heuristics - pattern match or nothing.
        
        NOTE: Previous version had overly broad fallback that caught ATP hydrolases
        and transport reactions. Precision improved by removing this.
        """
        # Must have water as reactant
        # Support multiple SMILES encodings: "O", "[H]O[H]", "[OH2]"
        WATER_SMILES = {"O", "[H]O[H]", "[OH2]", "O([H])[H]"}
        has_water = any(
            p.chebi_id == CHEBI_H2O or
            p.smiles in WATER_SMILES
            for p in reaction.left_participants
        )
        
        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant - not hydrolysis"
            )
        
        # ONLY pattern matching - no heuristic fallbacks
        # Use Hydrolase.PATTERNS explicitly (not self.PATTERNS) so subclasses
        # calling super() use the parent's generic patterns, not their own specific ones
        match = match_patterns(reaction, Hydrolase.PATTERNS)
        
        if match and match.matched:
            explanation = "Hydrolase: water-mediated bond cleavage"
            if match.bindings:
                substrate = match.bindings.get("substrate")
                if substrate:
                    explanation += f" of {substrate}"
            return ClassificationResult(
                is_member=True,
                explanation=explanation
            )
        
        # No pattern match = not hydrolase (removed overly broad fallback)
        return ClassificationResult(
            is_member=False,
            explanation="Pattern doesn't match expected hydrolysis transformation"
        )
