"""Oxidoreductase reducing C-O-C group reaction classification.

Oxidoreductases that reduce C-O-C (ether) bonds as acceptor (EC 1.23).
This includes lytic polysaccharide monooxygenases (LPMOs) that cleave
glycosidic bonds via an oxidative mechanism targeting C-O-C linkages.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseReducingCOCGroup(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which a reducing C-O-C group acts as acceptor."""

    GO_ID = "GO:0120546"  # oxidoreductase activity, reducing C-O-C group as acceptor
    EC_NUMBER_PREFIX = "1.23.-.-"  # Oxidoreductases reducing C-O-C group

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on C-O-C groups.

        Strategy:
        1. Must be an oxidoreductase (parent class)
        2. Look for ether bond (C-O-C) indicators in reaction labels
        3. Look for LPMO-related terms
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        # Look for ether bond or LPMO indicators in label
        ether_indicators = [
            "c-o-c",
            "ether",
            "lytic polysaccharide monooxygenase",
            "lpmo",
        ]

        for indicator in ether_indicators:
            if indicator in label_lower:
                return ClassificationResult(
                    is_member=True,
                    explanation=f"C-O-C oxidoreductase: label contains '{indicator}'",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No C-O-C / ether bond indicators detected",
        )
