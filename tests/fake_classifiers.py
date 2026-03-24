"""Fake classifier classes for testing.

These classifiers have simple, predictable logic that won't change,
allowing us to test the framework without being tied to real chemistry.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass


class FakeHydrolase(ReactionClass):
    """Test classifier that checks for 'FAKE_WATER' participant."""

    GO_ID = "GO:TEST001"
    EC_NUMBER_PREFIX = "3.99.99.99"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        # Simple check: has participant with smiles='FAKE_WATER'
        has_fake_water = any(
            p.smiles == "FAKE_WATER" for p in reaction.left_participants
        )

        if has_fake_water:
            return ClassificationResult(
                is_member=True, explanation="Has FAKE_WATER reactant"
            )
        return ClassificationResult(
            is_member=False, explanation="No FAKE_WATER reactant"
        )


class FakeOxidoreductase(ReactionClass):
    """Test classifier that checks for 'FAKE_NAD' participant."""

    GO_ID = "GO:TEST002"
    EC_NUMBER_PREFIX = "1.99.99.99"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        # Simple check: has participant with smiles containing 'FAKE_NAD'
        has_fake_nad = any(
            p.smiles and "FAKE_NAD" in p.smiles for p in reaction.all_participants()
        )

        if has_fake_nad:
            return ClassificationResult(
                is_member=True, explanation="Has FAKE_NAD cofactor"
            )
        return ClassificationResult(is_member=False, explanation="No FAKE_NAD cofactor")


class FakeTransferase(ReactionClass):
    """Test classifier that checks for equal participant count."""

    GO_ID = "GO:TEST003"
    EC_NUMBER_PREFIX = "2.99.99.99"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        # Simple check: same number of participants on each side
        if len(reaction.left_participants) == len(reaction.right_participants):
            return ClassificationResult(
                is_member=True, explanation="Equal participant count"
            )
        return ClassificationResult(
            is_member=False, explanation="Unequal participant count"
        )


class FakeLyase(ReactionClass):
    """Test classifier that checks for fragmentation pattern."""

    GO_ID = "GO:TEST004"
    EC_NUMBER_PREFIX = "4.99.99.99"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        # Simple check: fewer reactants than products
        if len(reaction.left_participants) < len(reaction.right_participants):
            return ClassificationResult(
                is_member=True, explanation="Fragmentation pattern"
            )
        return ClassificationResult(
            is_member=False, explanation="No fragmentation pattern"
        )


class AlwaysTrue(ReactionClass):
    """Test classifier that always returns True."""

    GO_ID = "GO:TEST999"
    EC_NUMBER_PREFIX = "99.99.99.99"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return ClassificationResult(
            is_member=True, explanation="Always true for testing"
        )


class AlwaysFalse(ReactionClass):
    """Test classifier that always returns False."""

    GO_ID = "GO:TEST000"
    EC_NUMBER_PREFIX = "0.0.0.0"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return ClassificationResult(
            is_member=False, explanation="Always false for testing"
        )
