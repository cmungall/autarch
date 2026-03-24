"""Simplified Kinase reaction classification using pattern DSL.

Kinases transfer phosphate groups from ATP/GTP to substrates.
Uses ultra-simple pattern matching.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import atp, adp, gtp, gdp, h_plus, p


class Kinase(ReactionClass):
    """kinase

    Core patterns:
    - ATP + substrate → ADP + phospho-substrate
    - GTP + substrate → GDP + phospho-substrate
    """

    GO_ID = "GO:0016301"  # kinase activity
    EC_NUMBER_PREFIX = "2.7.-.-"  # EC 2.7 = Transferases transferring phosphorus-containing groups (kinases)

    # Simple, declarative kinase patterns
    # ATP/GTP + substrate → ADP/GDP + phosphorylated product (+ optional H+)
    
    PATTERNS: list[Reaction] = [
        # ATP-based phosphorylation
        p(atp) + var("substrate") + optional(p(h_plus)) >> p(adp) + var("product") + optional(p(h_plus)),
        # GTP-based phosphorylation  
        p(gtp) + var("substrate") + optional(p(h_plus)) >> p(gdp) + var("product") + optional(p(h_plus)),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a kinase.

        Ultra-simple:
        1. Not a transport reaction
        2. Matches kinase pattern (ATP→ADP or GTP→GDP with substrate)
        """
        # Exclude transport reactions (ATPases that move molecules across membranes)
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction (ATPase/translocase) - not kinase"
            )

        # Pattern match for kinase reactions
        # Use strict matching - pattern should be closed (no extra participants)
        match = match_patterns(reaction, self.PATTERNS, strict=True)

        if match and match.matched:
            # Build explanation based on which pattern matched
            if "atp" in str(match.pattern).lower():
                explanation = "Kinase: ATP → ADP phosphoryl transfer"
            elif "gtp" in str(match.pattern).lower():
                explanation = "Kinase: GTP → GDP phosphoryl transfer"
            else:
                explanation = "Kinase: phosphoryl transfer"

            return ClassificationResult(
                is_member=True,
                explanation=explanation
            )

        return ClassificationResult(
            is_member=False,
            explanation="No kinase pattern found"
        )
