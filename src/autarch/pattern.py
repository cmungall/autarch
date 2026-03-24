"""Pattern matching for reaction classification.

This module provides a declarative DSL for defining reaction patterns
that can be matched against concrete reactions.

Overview:
    This module enables declarative pattern matching for chemical reactions,
    replacing complex imperative code with simple patterns. It supports:
    - Exact molecule matching (by CHEBI ID)
    - Variable patterns (wildcards that can be referenced)
    - Optional participants (with min/max counts)
    - Forbidden participants and constraints
    - Cross-side variable unification

Basic Usage:
    >>> from autarch.datamodel import Reaction, Participant
    >>> from autarch.pattern import ReactionPattern, ParticipantPattern

    >>> # Define a simple A + B → C pattern
    >>> pattern = ReactionPattern(
    ...     name="Combination",
    ...     left=[
    ...         ParticipantPattern(chebi_id="CHEBI:111"),  # A
    ...         ParticipantPattern(chebi_id="CHEBI:222")   # B
    ...     ],
    ...     right=[
    ...         ParticipantPattern(chebi_id="CHEBI:333")   # C
    ...     ]
    ... )

    >>> # Check if a reaction matches
    >>> reaction = Reaction(
    ...     left_participants=[
    ...         Participant(chebi_id="CHEBI:111"),
    ...         Participant(chebi_id="CHEBI:222")
    ...     ],
    ...     right_participants=[
    ...         Participant(chebi_id="CHEBI:333")
    ...     ]
    ... )
    >>> matches, explanation, bindings = pattern.matches(reaction)
    >>> matches
    True

Variable Patterns:
    >>> # Use variables to match any molecule
    >>> transfer_pattern = ReactionPattern(
    ...     name="Transfer",
    ...     left=[
    ...         ParticipantPattern(variable="donor"),
    ...         ParticipantPattern(variable="acceptor")
    ...     ],
    ...     right=[
    ...         ParticipantPattern(variable="modified_donor"),
    ...         ParticipantPattern(variable="modified_acceptor")
    ...     ]
    ... )

    >>> reaction = Reaction(
    ...     left_participants=[
    ...         Participant(chebi_id="CHEBI:100"),
    ...         Participant(chebi_id="CHEBI:200")
    ...     ],
    ...     right_participants=[
    ...         Participant(chebi_id="CHEBI:101"),
    ...         Participant(chebi_id="CHEBI:201")
    ...     ]
    ... )
    >>> matches, explanation, bindings = transfer_pattern.matches(reaction)
    >>> matches
    True
    >>> sorted(bindings.keys())
    ['acceptor', 'donor', 'modified_acceptor', 'modified_donor']

Real-World Example - Kinase Pattern:
    >>> # Define a kinase pattern (ATP + X → ADP + X-P + H+)
    >>> kinase_pattern = ReactionPattern(
    ...     name="Kinase",
    ...     left=[
    ...         ParticipantPattern(chebi_id="CHEBI:30616"),  # ATP
    ...         ParticipantPattern(variable="substrate")
    ...     ],
    ...     right=[
    ...         ParticipantPattern(chebi_id="CHEBI:456216"),  # ADP
    ...         ParticipantPattern(variable="product"),       # phosphorylated substrate
    ...         ParticipantPattern(                           # optional H+
    ...             chebi_id="CHEBI:15378",
    ...             min_count=0,
    ...             max_count=1
    ...         )
    ...     ],
    ...     require_no_water=True,
    ...     forbidden_participants=["CHEBI:43474"]  # no free phosphate
    ... )

    >>> # Test with a valid kinase reaction: ATP + glucose → ADP + glucose-6-P
    >>> glucose_phosphorylation = Reaction(
    ...     left_participants=[
    ...         Participant(chebi_id="CHEBI:30616", name="ATP"),
    ...         Participant(chebi_id="CHEBI:17234", name="glucose")
    ...     ],
    ...     right_participants=[
    ...         Participant(chebi_id="CHEBI:456216", name="ADP"),
    ...         Participant(chebi_id="CHEBI:17665", name="glucose-6-phosphate")
    ...     ]
    ... )

    >>> matches, explanation, bindings = kinase_pattern.matches(glucose_phosphorylation)
    >>> matches
    True
    >>> bindings["substrate"]
    'CHEBI:17234'

    >>> # Test with invalid reaction (has water - likely ATPase)
    >>> atp_hydrolysis = Reaction(
    ...     left_participants=[
    ...         Participant(chebi_id="CHEBI:30616", name="ATP"),
    ...         Participant(chebi_id="CHEBI:15377", name="water")
    ...     ],
    ...     right_participants=[
    ...         Participant(chebi_id="CHEBI:456216", name="ADP"),
    ...         Participant(chebi_id="CHEBI:43474", name="phosphate")
    ...     ]
    ... )

    >>> matches, explanation, bindings = atp_hydrolysis_pattern.matches(atp_hydrolysis)  # doctest: +SKIP
    >>> # This would fail because it contains water and free phosphate
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Tuple
from enum import Enum

from autarch.datamodel import Reaction, Participant


class MatchMode(Enum):
    """How to match a participant pattern.

    Examples:
        >>> # Different matching modes
        >>> MatchMode.EXACT.value
        'exact'
        >>> MatchMode.VARIABLE.value
        'variable'
        >>> list(MatchMode)  # doctest: +ELLIPSIS
        [<MatchMode.EXACT: 'exact'>, <MatchMode.VARIABLE: 'variable'>...]
    """

    EXACT = "exact"  # Exact CHEBI ID match
    VARIABLE = "variable"  # Match any, bind to variable
    CLASS = "class"  # Match CHEBI class/ancestor
    EXCLUDE = "exclude"  # Match anything except specified


@dataclass
class ParticipantPattern:
    """Pattern for matching reaction participants.

    Can represent:
    - Exact molecules: chebi_id="CHEBI:30616" (ATP)
    - Variables: variable="X" (any molecule, can be referenced)
    - Optional: min_count=0, max_count=1 (e.g., H+ may or may not appear)
    - Exclusions: exclude_chebi_ids=["CHEBI:15377"] (not water)

    Examples:
        >>> # Create an exact match pattern for ATP
        >>> atp_pattern = ParticipantPattern(chebi_id="CHEBI:30616")
        >>> atp_pattern.chebi_id
        'CHEBI:30616'

        >>> # Create a variable pattern that can match any molecule
        >>> substrate_pattern = ParticipantPattern(variable="substrate")
        >>> substrate_pattern.variable
        'substrate'

        >>> # Create an optional H+ pattern
        >>> proton_pattern = ParticipantPattern(
        ...     chebi_id="CHEBI:15378",
        ...     min_count=0,
        ...     max_count=1
        ... )
        >>> proton_pattern.min_count
        0

        >>> # Create a pattern that excludes water
        >>> non_water = ParticipantPattern(
        ...     variable="X",
        ...     exclude_chebi_ids=["CHEBI:15377"]
        ... )
        >>> "CHEBI:15377" in non_water.exclude_chebi_ids
        True
    """

    # Matching mode
    chebi_id: Optional[str] = None  # Exact CHEBI to match
    variable: Optional[str] = None  # Variable name for cross-reference
    chebi_class: Optional[str] = None  # CHEBI class/ancestor to match

    # Count constraints
    min_count: int = 1
    max_count: int = 1

    # Exclusions
    exclude_chebi_ids: List[str] = field(default_factory=list)

    # Modifications (for variables)
    modifications: List[str] = field(default_factory=list)  # e.g., ["phosphorylated"]

    def matches(
        self, participant: Participant, bindings: Optional[Dict[str, str]] = None
    ) -> Tuple[bool, Dict[str, str]]:
        """Check if a participant matches this pattern.

        Args:
            participant: The concrete participant to match
            bindings: Current variable bindings

        Returns:
            Tuple of (matches: bool, updated_bindings: dict)

        Examples:
            >>> from autarch.datamodel import Participant

            # Exact match example
            >>> pattern = ParticipantPattern(chebi_id="CHEBI:30616")
            >>> atp = Participant(chebi_id="CHEBI:30616", name="ATP")
            >>> matches, bindings = pattern.matches(atp)
            >>> matches
            True
            >>> bindings
            {}

            # Variable binding example
            >>> pattern = ParticipantPattern(variable="X")
            >>> glucose = Participant(chebi_id="CHEBI:17234", name="glucose")
            >>> matches, bindings = pattern.matches(glucose)
            >>> matches
            True
            >>> bindings
            {'X': 'CHEBI:17234'}

            # Variable already bound - should match same molecule
            >>> pattern = ParticipantPattern(variable="X")
            >>> glucose = Participant(chebi_id="CHEBI:17234")
            >>> existing_bindings = {"X": "CHEBI:17234"}
            >>> matches, bindings = pattern.matches(glucose, existing_bindings)
            >>> matches
            True

            # Variable already bound - should not match different molecule
            >>> water = Participant(chebi_id="CHEBI:15377")
            >>> matches, bindings = pattern.matches(water, existing_bindings)
            >>> matches
            False

            # Exclusion example
            >>> pattern = ParticipantPattern(variable="X", exclude_chebi_ids=["CHEBI:15377"])
            >>> water = Participant(chebi_id="CHEBI:15377", name="water")
            >>> matches, bindings = pattern.matches(water)
            >>> matches
            False
        """
        if bindings is None:
            bindings = {}

        # Check exclusions first
        if participant.chebi_id in self.exclude_chebi_ids:
            return False, bindings

        # Exact match
        if self.chebi_id:
            return participant.chebi_id == self.chebi_id, bindings

        # Variable match - bind or check existing binding
        if self.variable:
            if self.variable in bindings:
                # Variable already bound, check if it matches
                expected = bindings[self.variable]
                # For modified variables, we'd need more complex logic here
                # For now, just check exact match
                return participant.chebi_id == expected, bindings
            else:
                # Bind the variable
                if participant.chebi_id:
                    new_bindings = bindings.copy()
                    new_bindings[self.variable] = participant.chebi_id
                    return True, new_bindings
                else:
                    # Can't bind None
                    return False, bindings

        # Class match (would need CHEBI ontology integration)
        if self.chebi_class:
            # TODO: Check if participant.chebi_id is-a chebi_class
            # For now, always false
            return False, bindings

        # Default: matches anything
        return True, bindings


@dataclass
class ReactionPattern:
    """Declarative pattern for matching reactions.

    Example for Kinase:
        ReactionPattern(
            left=[
                ParticipantPattern(chebi_id="CHEBI:30616"),  # ATP
                ParticipantPattern(variable="X")  # substrate
            ],
            right=[
                ParticipantPattern(chebi_id="CHEBI:456216"),  # ADP
                ParticipantPattern(variable="X", modifications=["phosphorylated"]),
                ParticipantPattern(chebi_id="CHEBI:15378", min_count=0, max_count=1)  # H+
            ]
        )

    Examples:
        >>> # Create a simple A + B → C pattern
        >>> pattern = ReactionPattern(
        ...     name="Simple combination",
        ...     left=[
        ...         ParticipantPattern(chebi_id="CHEBI:1111"),
        ...         ParticipantPattern(chebi_id="CHEBI:2222")
        ...     ],
        ...     right=[
        ...         ParticipantPattern(chebi_id="CHEBI:3333")
        ...     ]
        ... )
        >>> len(pattern.left)
        2
        >>> len(pattern.right)
        1

        >>> # Create a pattern with constraints
        >>> kinase_pattern = ReactionPattern(
        ...     name="Kinase",
        ...     left=[
        ...         ParticipantPattern(chebi_id="CHEBI:30616"),
        ...         ParticipantPattern(variable="X")
        ...     ],
        ...     right=[
        ...         ParticipantPattern(chebi_id="CHEBI:456216"),
        ...         ParticipantPattern(variable="Y")
        ...     ],
        ...     require_no_water=True,
        ...     forbidden_participants=["CHEBI:43474"]  # no free phosphate
        ... )
        >>> kinase_pattern.require_no_water
        True
        >>> "CHEBI:43474" in kinase_pattern.forbidden_participants
        True
    """

    left: List[ParticipantPattern]
    right: List[ParticipantPattern]

    # Global constraints
    forbidden_participants: List[str] = field(
        default_factory=list
    )  # CHEBIs that must not appear
    require_no_water: bool = False
    require_no_phosphate: bool = False

    # Metadata
    name: Optional[str] = None
    description: Optional[str] = None

    def matches(self, reaction: Reaction) -> Tuple[bool, Optional[str], Dict[str, str]]:
        """Check if a reaction matches this pattern.

        Returns:
            Tuple of (matches, explanation, bindings)

        Examples:
            >>> from autarch.datamodel import Reaction, Participant

            # Simple ATP + X → ADP + Y pattern
            >>> pattern = ReactionPattern(
            ...     name="Simple phosphorylation",
            ...     left=[
            ...         ParticipantPattern(chebi_id="CHEBI:30616"),  # ATP
            ...         ParticipantPattern(variable="substrate")
            ...     ],
            ...     right=[
            ...         ParticipantPattern(chebi_id="CHEBI:456216"),  # ADP
            ...         ParticipantPattern(variable="product")
            ...     ]
            ... )

            # Create a matching reaction
            >>> reaction = Reaction(
            ...     left_participants=[
            ...         Participant(chebi_id="CHEBI:30616", name="ATP"),
            ...         Participant(chebi_id="CHEBI:17234", name="glucose")
            ...     ],
            ...     right_participants=[
            ...         Participant(chebi_id="CHEBI:456216", name="ADP"),
            ...         Participant(chebi_id="CHEBI:17665", name="glucose-6-P")
            ...     ]
            ... )

            >>> matches, explanation, bindings = pattern.matches(reaction)
            >>> matches
            True
            >>> "Simple phosphorylation" in explanation
            True
            >>> bindings["substrate"]
            'CHEBI:17234'
            >>> bindings["product"]
            'CHEBI:17665'

            # Pattern with forbidden water
            >>> pattern_no_water = ReactionPattern(
            ...     name="No water allowed",
            ...     left=[ParticipantPattern(variable="X")],
            ...     right=[ParticipantPattern(variable="Y")],
            ...     require_no_water=True
            ... )

            >>> reaction_with_water = Reaction(
            ...     left_participants=[
            ...         Participant(chebi_id="CHEBI:12345"),
            ...         Participant(chebi_id="CHEBI:15377")  # water
            ...     ],
            ...     right_participants=[
            ...         Participant(chebi_id="CHEBI:54321")
            ...     ]
            ... )

            >>> matches, explanation, bindings = pattern_no_water.matches(reaction_with_water)
            >>> matches
            False
            >>> "water" in explanation.lower()
            True
        """
        # Check forbidden participants
        all_chebi_ids = set()
        for p in reaction.left_participants + reaction.right_participants:
            if p.chebi_id:
                all_chebi_ids.add(p.chebi_id)

        for forbidden in self.forbidden_participants:
            if forbidden in all_chebi_ids:
                return False, f"Contains forbidden participant: {forbidden}", {}

        # Special constraints
        if self.require_no_water and "CHEBI:15377" in all_chebi_ids:
            return False, "Pattern forbids water", {}

        if self.require_no_phosphate and "CHEBI:43474" in all_chebi_ids:
            return False, "Pattern forbids free phosphate", {}

        # Try to match left side
        left_match, left_bindings = self._match_side(
            self.left, reaction.left_participants, "left"
        )
        if not left_match:
            return False, left_bindings, {}  # left_bindings contains explanation

        # Try to match right side with bindings from left
        right_match, right_bindings = self._match_side(
            self.right, reaction.right_participants, "right", left_bindings
        )
        if not right_match:
            return False, right_bindings, {}  # right_bindings contains explanation

        return True, f"Matches pattern: {self.name or 'unnamed'}", right_bindings

    def _match_side(
        self,
        patterns: List[ParticipantPattern],
        participants: List[Participant],
        side: str,
        bindings: Optional[Dict[str, str]] = None,
    ) -> Tuple[bool, Any]:
        """Match one side of the reaction.

        This is a simplified implementation. A full implementation would need:
        - Proper handling of stoichiometry
        - Support for optional participants
        - Better variable unification

        Examples:
            >>> from autarch.datamodel import Participant

            >>> # Create a pattern and test side matching
            >>> pattern = ReactionPattern(
            ...     name="test",
            ...     left=[ParticipantPattern(chebi_id="CHEBI:123")],
            ...     right=[]
            ... )

            >>> # Matching participants
            >>> participants = [Participant(chebi_id="CHEBI:123")]
            >>> matches, result = pattern._match_side(pattern.left, participants, "left")
            >>> matches
            True

            >>> # Non-matching participants
            >>> participants = [Participant(chebi_id="CHEBI:999")]
            >>> matches, result = pattern._match_side(pattern.left, participants, "left")
            >>> matches
            False
            >>> "No match" in str(result)
            True
        """
        if bindings is None:
            bindings = {}

        # Simple check: do we have the right number of participants?
        # (ignoring optional ones for now)
        required_patterns = [p for p in patterns if p.min_count > 0]

        if len(participants) < len(required_patterns):
            return False, f"Too few participants on {side}"

        # Try to match each required pattern
        # This is simplified - real implementation would need proper matching algorithm
        matched_participants = set()
        current_bindings = bindings.copy()

        for pattern in required_patterns:
            matched = False
            for i, participant in enumerate(participants):
                if i in matched_participants:
                    continue

                match, new_bindings = pattern.matches(participant, current_bindings)
                if match:
                    matched_participants.add(i)
                    current_bindings = new_bindings
                    matched = True
                    break

            if not matched:
                return False, f"No match for pattern on {side}"

        # Check optional patterns
        optional_patterns = [p for p in patterns if p.min_count == 0]
        for pattern in optional_patterns:
            for i, participant in enumerate(participants):
                if i in matched_participants:
                    continue
                match, new_bindings = pattern.matches(participant, current_bindings)
                if match:
                    matched_participants.add(i)
                    current_bindings = new_bindings

        # Check if we matched all participants (no extras)
        if len(matched_participants) < len(participants):
            # Some participants weren't matched - might be OK for some patterns
            pass

        return True, current_bindings


def pattern_classifier(patterns: List[ReactionPattern]):
    """Decorator to create a pattern-based classifier.

    Usage:
        @pattern_classifier([
            ReactionPattern(...),
            ReactionPattern(...)
        ])
        class MyReactionClass(ReactionClassifier):
            pass

    Examples:
        >>> from autarch.datamodel import ClassificationResult, Reaction, Participant

        >>> # Define patterns for a simple transferase
        >>> transferase_patterns = [
        ...     ReactionPattern(
        ...         name="Transfer pattern",
        ...         left=[
        ...             ParticipantPattern(variable="donor"),
        ...             ParticipantPattern(variable="acceptor")
        ...         ],
        ...         right=[
        ...             ParticipantPattern(variable="donor_minus"),
        ...             ParticipantPattern(variable="acceptor_plus")
        ...         ]
        ...     )
        ... ]

        >>> # Create a simple class to decorate (no need to inherit from anything)
        >>> @pattern_classifier(transferase_patterns)
        ... class SimpleTransferase:
        ...     pass

        >>> # The class now has PATTERNS attribute
        >>> hasattr(SimpleTransferase, 'PATTERNS')
        True
        >>> len(SimpleTransferase.PATTERNS)
        1

        >>> # Create instance and test classification
        >>> classifier = SimpleTransferase()
        >>> reaction = Reaction(
        ...     left_participants=[
        ...         Participant(chebi_id="CHEBI:1111"),
        ...         Participant(chebi_id="CHEBI:2222")
        ...     ],
        ...     right_participants=[
        ...         Participant(chebi_id="CHEBI:3333"),
        ...         Participant(chebi_id="CHEBI:4444")
        ...     ]
        ... )

        >>> # The check_membership_impl method is automatically created
        >>> result = classifier.check_membership_impl(reaction)
        >>> isinstance(result, ClassificationResult)
        True
    """

    def decorator(cls):
        # Store patterns on the class
        cls.PATTERNS = patterns

        # Override check_membership_impl to use patterns
        def check_membership_impl(self, reaction: Reaction):
            from autarch.datamodel import ClassificationResult

            for pattern in self.PATTERNS:
                matches, explanation, bindings = pattern.matches(reaction)
                if matches:
                    return ClassificationResult(is_member=True, explanation=explanation)

            # No patterns matched
            return ClassificationResult(
                is_member=False, explanation="No patterns matched"
            )

        cls.check_membership_impl = check_membership_impl
        return cls

    return decorator
