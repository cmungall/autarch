"""Ligase forming nitrogen-metal bonds classification.

EC 6.6 ligases catalyze the formation of coordination bonds between
nitrogen atoms and metal ions, typically inserting metals into
tetrapyrrole rings (porphyrins, corrins).
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ligase import Ligase


class LigaseFormingNitrogenMetalBonds(Ligase):
    """ligase activity, forming nitrogen-metal bonds

    Catalysis of the joining of a metal ion and a nitrogen atom,
    with the concomitant hydrolysis of ATP or a similar triphosphate.

    Examples:
    - Magnesium chelatase: inserts Mg2+ into protoporphyrin IX (chlorophyll biosynthesis)
    - Ferrochelatase: inserts Fe2+ into protoporphyrin IX (heme biosynthesis)
    - Cobalt chelatase: inserts Co2+ into precorrin-2 (vitamin B12 biosynthesis)
    - Nickel chelatase: inserts Ni2+ into sirohydrochlorin (coenzyme F430)
    """

    GO_ID = "GO:0051002"  # ligase activity, forming nitrogen-metal bonds
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "6.6.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a nitrogen-metal bond-forming ligase.

        Strategy:
        1. Must be a ligase (ATP-dependent bond formation)
        2. Look for chelatase-related terms by label
        3. Look for porphyrin/tetrapyrrole substrates
        """
        # First check if it's a ligase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a ligase: {parent_result.explanation}",
            )

        # Label-based detection for chelatases
        all_names = [
            p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
            if p.name
        ]
        chelatase_terms = {
            "chelatase", "magnesium chelatase", "cobalt chelatase",
            "ferrochelatase", "nickel chelatase",
        }
        if any(term in name for name in all_names for term in chelatase_terms):
            return ClassificationResult(
                is_member=True,
                explanation="N-metal ligase: chelatase detected by label",
            )

        # Check for porphyrin/tetrapyrrole substrates
        tetrapyrrole_terms = {
            "porphyrin", "protoporphyrin", "corrin", "precorrin",
            "sirohydrochlorin", "tetrapyrrole",
        }
        if any(term in name for name in all_names for term in tetrapyrrole_terms):
            return ClassificationResult(
                is_member=True,
                explanation="N-metal ligase: tetrapyrrole substrate with ATP-dependent metal insertion",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No chelatase/tetrapyrrole pattern - not a nitrogen-metal bond-forming ligase",
        )
