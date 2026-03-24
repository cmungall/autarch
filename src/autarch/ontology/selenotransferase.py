"""Selenotransferase reaction classification.

Selenotransferases (EC 2.9) transfer selenium-containing groups
from one compound to another. Key substrates include selenophosphate
and selenocysteine.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase


class Selenotransferase(Transferase):
    """Catalysis of the transfer of a selenium-containing group from one compound (donor) to another (acceptor)."""

    GO_ID = "GO:0016785"  # selenotransferase activity
    EC_NUMBER_PREFIX = "2.9.-.-"  # Transferring selenium-containing groups

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a selenotransferase.

        Strategy:
        1. Must be a transferase (parent class) OR have selenium indicators
        2. Look for selenium-containing participants (Se in SMILES)
        3. Look for selenium-related terms in labels
        """
        label_lower = reaction.label.lower() if reaction.label else ""

        # Look for selenium indicators in label
        selenium_label_indicators = [
            "seleno",
            "selenium",
            "selenide",
        ]

        has_selenium_label = any(
            indicator in label_lower for indicator in selenium_label_indicators
        )

        # Look for selenium in participant SMILES
        has_selenium_smiles = any(
            p.smiles and "[Se" in p.smiles
            for p in reaction.left_participants + reaction.right_participants
        )

        if not has_selenium_label and not has_selenium_smiles:
            return ClassificationResult(
                is_member=False,
                explanation="No selenium-containing groups detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Selenotransferase: selenium-containing group transfer",
        )
