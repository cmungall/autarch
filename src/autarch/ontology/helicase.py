"""Helicase activity classification.

Helicases unwind DNA/RNA double helices using ATP hydrolysis.
GO molecular function, EC 3.6.4.x classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.pattern_dsl import var, optional
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_H2O, CHEBI_PHOSPHATE, adp, atp, phosphate, h_plus, p


class Helicase(ReactionClass):
    """helicase

    Examples:
    - DNA helicase: dsDNA + ATP + H2O → 2 ssDNA + ADP + Pi
    - RNA helicase: dsRNA + ATP + H2O → 2 ssRNA + ADP + Pi
    - RecA helicase: DNA duplex + ATP → unwound DNA + ADP + Pi
    - Replicative helicases during DNA replication
    """

    GO_ID = "GO:0004386"  # helicase activity
    EC_NUMBER_PREFIX = "3.6.4.-"  # Acting on acid anhydrides; involved in cellular processes

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # DNA helicase: dsDNA + ATP + H2O → ssDNA + ssDNA + ADP + Pi
        var("double_stranded") + p(atp) + var("water") >> var("strand1") + var("strand2") + p(adp) + p(phosphate) + optional(h_plus),
        # RNA helicase: dsRNA + ATP → ssRNA + ssRNA + ADP + Pi
        var("dsrna") + p(atp) >> var("ssrna1") + var("ssrna2") + p(adp) + p(phosphate),
        # General nucleic acid unwinding: NA duplex + ATP → separated strands + ADP + Pi
        var("duplex") + p(atp) + var("cofactor") >> var("separated") + p(adp) + p(phosphate) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction involves helicase activity.

        Strategy:
        1. Must be ATP hydrolysis activity (parent class)
        2. Must involve nucleic acid substrates (DNA/RNA)
        3. Must involve strand separation/unwinding
        4. Look for helicase-specific patterns
        """
        left_ids = {participant.chebi_id for participant in reaction.left_participants}
        right_ids = {participant.chebi_id for participant in reaction.right_participants}
        has_atp_hydrolysis = (
            CHEBI_ATP in left_ids
            and CHEBI_H2O in left_ids
            and CHEBI_ADP in right_ids
            and CHEBI_PHOSPHATE in right_ids
        )
        if not has_atp_hydrolysis:
            return ClassificationResult(
                is_member=False,
                explanation="Missing ATP hydrolysis signature required for helicase activity",
            )

        # Helicase requires nucleic acid substrates - currently no ChEBI-based detection
        # This class will not match any reactions until nucleic acid detection is implemented
        return ClassificationResult(
            is_member=False,
            explanation="Helicase detection requires nucleic acid substrate identification"
        )
