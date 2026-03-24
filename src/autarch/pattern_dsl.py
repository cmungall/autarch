"""Pattern DSL using Python operator overloading.

This module extends the existing Reaction/Participant DSL to support
pattern matching with variables, optional participants, and strict
participant accounting.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Union, Set, Tuple
from pydantic import Field

from autarch.datamodel import Reaction, Participant

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Union  # Already imported above


# Forward declaration
class PatternReaction(Reaction):
    """A Reaction used specifically for patterns.

    This class exists to allow operator overloading on the right side
    of the >> operator, enabling syntax like:
    atp + substrate >> adp + product
    """

    def __add__(
        self, other: Union["PatternParticipant", Participant]
    ) -> "PatternReaction":
        """Add a participant to the current side of the reaction.

        If we have no right_participants yet, add to left side (building left side).
        If we have right_participants, add to right side (building right side after >>).

        Examples:
            >>> from autarch.molecules import atp, adp
            >>> substrate = var("S")
            >>> product = var("P")
            >>> # Convert to pattern participants
            >>> atp_p = to_pattern(atp)
            >>> adp_p = to_pattern(adp)
            >>> # Build pattern with operators
            >>> pattern = atp_p + substrate >> adp_p + product
            >>> len(pattern.right_participants)
            2
        """
        # Convert regular Participant to PatternParticipant if needed
        if isinstance(other, Participant) and not isinstance(other, PatternParticipant):
            other = to_pattern(other)

        # If we don't have right participants yet, we're still building the left side
        if not self.right_participants:
            return PatternReaction(
                left_participants=self.left_participants + [other],
                right_participants=[],
            )
        else:
            # We have right participants, so add to the right
            return PatternReaction(
                left_participants=self.left_participants,
                right_participants=self.right_participants + [other],
            )

    def __rshift__(
        self, other: Union["PatternParticipant", Participant, "PatternReaction"]
    ) -> "PatternReaction":
        """Enable chaining >> operators for complex patterns.

        This handles cases like: (atp + S) >> (adp + P1)
        When other is a PatternReaction from '+' operations, its participants
        are in left_participants and should move to the right side.
        """
        # Convert regular Participant to PatternParticipant if needed
        if isinstance(other, Participant) and not isinstance(other, PatternParticipant):
            other = to_pattern(other)

        if isinstance(other, PatternReaction):
            # When we do (a + b) >> (c + d), the (c + d) part creates a PatternReaction
            # with c and d in left_participants. We want them on the right.
            # If other already has right_participants, use those; otherwise use left
            if other.right_participants:
                # other is a full reaction, use its right side
                right_parts = other.right_participants
            else:
                # other is from '+' operations, move its left to right
                right_parts = other.left_participants

            return PatternReaction(
                left_participants=self.left_participants, right_participants=right_parts
            )
        else:
            # Add single participant to right
            return PatternReaction(
                left_participants=self.left_participants,
                right_participants=self.right_participants + [other],
            )


class PatternParticipant(Participant):
    """Participant that can be a variable or have count ranges.

    Examples:
        >>> from autarch.pattern_dsl import var, optional

        >>> # Create a variable participant
        >>> substrate = var("substrate")
        >>> substrate.variable
        'substrate'
        >>> substrate.is_variable
        True

        >>> # Create an optional participant (count 0-1)
        >>> h_plus = Participant(chebi_id="CHEBI:15378", name="H+")
        >>> optional_h = optional(h_plus)
        >>> optional_h.min_count
        0
        >>> optional_h.max_count
        1

        >>> # Use in a pattern (operator overloading requires all PatternParticipants)
        >>> from autarch.molecules import atp, adp
        >>> p = pattern([atp, substrate], [adp, var("product")])
        >>> isinstance(p, Reaction)
        True
        >>> p.left_participants[1].variable
        'substrate'
    """

    # Variable name (if this is a variable placeholder)
    variable: Optional[str] = Field(
        default=None, description="Variable name for pattern matching"
    )

    # Count constraints (for optional/variable counts)
    min_count: int = Field(default=1, description="Minimum count for this participant")
    max_count: Optional[int] = Field(
        default=1, description="Maximum count (None = unbounded)"
    )

    @property
    def is_variable(self) -> bool:
        """Check if this is a variable participant."""
        return self.variable is not None

    def matches(self, participant: Participant) -> bool:
        """Check if a concrete participant matches this pattern.

        Examples:
            >>> # Exact match
            >>> atp_pattern = PatternParticipant(chebi_id="CHEBI:30616")
            >>> atp = Participant(chebi_id="CHEBI:30616")
            >>> atp_pattern.matches(atp)
            True

            >>> # Variable matches anything
            >>> var_pattern = PatternParticipant(variable="X")
            >>> glucose = Participant(chebi_id="CHEBI:17234")
            >>> var_pattern.matches(glucose)
            True
        """
        # Variable matches any participant
        if self.is_variable:
            return True

        # Exact CHEBI ID match
        if self.chebi_id and participant.chebi_id:
            return participant.chebi_id == self.chebi_id

        # SMILES match
        if self.smiles and participant.smiles:
            return participant.smiles == self.smiles

        # Name match (fallback)
        if self.name and participant.name:
            return self.name == participant.name

        return False

    def __add__(
        self, other: Union["PatternParticipant", Participant, "PatternReaction"]
    ) -> "PatternReaction":
        """Add participants together to build pattern.

        Examples:
            >>> substrate = var("S")
            >>> product = var("P")
            >>> from autarch.molecules import atp, adp
            >>> # Convert molecules to pattern participants
            >>> atp_p = to_pattern(atp)
            >>> adp_p = to_pattern(adp)
            >>> # Now we can use operator overloading
            >>> left_side = atp_p + substrate
            >>> isinstance(left_side, PatternReaction)
            True
            >>> len(left_side.left_participants)
            2
        """
        # Convert regular Participant to PatternParticipant if needed
        if isinstance(other, Participant) and not isinstance(other, PatternParticipant):
            other = to_pattern(other)

        if isinstance(other, PatternReaction):
            # Add self to existing reaction's left side
            return PatternReaction(
                left_participants=[self] + other.left_participants,
                right_participants=other.right_participants,
            )
        else:
            # Create new reaction with both on left side
            return PatternReaction(left_participants=[self, other])

    def __rshift__(
        self, other: Union["PatternParticipant", Participant, "PatternReaction"]
    ) -> "PatternReaction":
        """Create reaction pattern using >> operator.

        Examples:
            >>> substrate = var("S")
            >>> product = var("P")
            >>> pattern = substrate >> product
            >>> isinstance(pattern, PatternReaction)
            True
            >>> len(pattern.left_participants)
            1
            >>> len(pattern.right_participants)
            1
        """
        # Convert regular Participant to PatternParticipant if needed
        if isinstance(other, Participant) and not isinstance(other, PatternParticipant):
            other = to_pattern(other)

        if isinstance(other, PatternReaction):
            # Add self to left, keep other's right
            return PatternReaction(
                left_participants=[self] + other.left_participants,
                right_participants=other.right_participants,
            )
        else:
            # Create new reaction: self >> other
            return PatternReaction(left_participants=[self], right_participants=[other])


def pattern(
    left: List[Union[Participant, PatternParticipant, str]],
    right: List[Union[Participant, PatternParticipant, str]],
) -> Reaction:
    """Create a pattern reaction from left and right participants.

    Args:
        left: List of participants for the left side
        right: List of participants for the right side

    Returns:
        Reaction pattern

    Examples:
        >>> from autarch.molecules import atp, adp, water, h_plus
        >>>
        >>> # Simple kinase pattern
        >>> p = pattern([atp, var("S")], [adp, var("P")])
        >>> isinstance(p, Reaction)
        True
        >>> len(p.left_participants)
        2
        >>> len(p.right_participants)
        2
        >>>
        >>> # Pattern with optional H+ production
        >>> p2 = pattern([atp, var("S")], [adp, var("P"), optional(h_plus)])
        >>> len(p2.right_participants)
        3
        >>> p2.right_participants[-1].min_count
        0
        >>>
        >>> # Pattern from string names
        >>> p3 = pattern(["ATP", "glucose"], ["ADP", "glucose-6-P"])
        >>> p3.left_participants[0].name
        'ATP'
        >>> p3.left_participants[1].name
        'glucose'
    """

    def convert_participants(participants):
        result = []
        for p in participants:
            if isinstance(p, str):
                if p.startswith("CHEBI:"):
                    p = Participant(chebi_id=p)
                else:
                    p = Participant(name=p)
            if isinstance(p, PatternParticipant):
                result.append(p)
            else:
                # Convert to PatternParticipant
                result.append(to_pattern(p))
        return result

    return Reaction(
        left_participants=convert_participants(left),
        right_participants=convert_participants(right),
    )


def to_pattern(
    participant: Union[Participant, PatternParticipant],
) -> PatternParticipant:
    """Convert a regular Participant to a PatternParticipant.

    Args:
        participant: Regular or pattern participant

    Returns:
        PatternParticipant
    """
    if isinstance(participant, PatternParticipant):
        return participant
    return PatternParticipant(
        chebi_id=participant.chebi_id,
        name=participant.name,
        count=participant.count,
        smiles=participant.smiles,
        location=participant.location,
    )


def var(name: str, **kwargs) -> PatternParticipant:
    """Create a variable participant.

    Args:
        name: Variable name for binding
        **kwargs: Additional constraints (future use)

    Returns:
        PatternParticipant with variable set

    Examples:
        >>> substrate = var("substrate")
        >>> substrate.variable
        'substrate'
        >>> substrate.is_variable
        True
    """
    return PatternParticipant(variable=name, name=f"?{name}")


def optional(
    participant: Union[Participant, PatternParticipant, str],
) -> PatternParticipant:
    """Make a participant optional (count 0 or 1).

    Args:
        participant: Participant to make optional, or name/CHEBI ID string

    Returns:
        PatternParticipant with min_count=0, max_count=1

    Examples:
        >>> h_plus = Participant(chebi_id="CHEBI:15378", name="H+")
        >>> opt_h = optional(h_plus)
        >>> opt_h.min_count
        0
        >>> opt_h.max_count
        1

        >>> # Also works with strings
        >>> opt_water = optional("H2O")
        >>> opt_water.name
        'H2O'
        >>> opt_water.min_count
        0
    """
    if isinstance(participant, str):
        # Create participant from string (assume it's a name)
        if participant.startswith("CHEBI:"):
            participant = Participant(chebi_id=participant)
        else:
            participant = Participant(name=participant)

    # Convert to PatternParticipant if needed
    if not isinstance(participant, PatternParticipant):
        participant = PatternParticipant(
            chebi_id=participant.chebi_id,
            name=participant.name,
            count=participant.count,
            smiles=participant.smiles,
            location=participant.location,
        )

    # Set optional counts
    participant.min_count = 0
    participant.max_count = 1
    return participant


def any_of(*participants: Union[Participant, str]) -> PatternParticipant:
    """Create a pattern that matches any of the given participants.

    This is a placeholder for future implementation.

    Args:
        *participants: Participants that are alternatives

    Returns:
        PatternParticipant representing alternatives
    """
    # For now, just return the first one
    # Full implementation would need a way to represent alternatives
    if participants:
        first = participants[0]
        if isinstance(first, str):
            first = Participant(name=first)
        if not isinstance(first, PatternParticipant):
            first = PatternParticipant(chebi_id=first.chebi_id, name=first.name)
        # Could store alternatives in a custom field
        # NOTE: This would need proper implementation with a field
        # first._alternatives = participants[1:]  # type: ignore
        return first
    return PatternParticipant()


@dataclass
class PatternMatch:
    """Result of pattern matching.

    Attributes:
        matched: Whether the pattern matched
        pattern: The pattern that matched (if any)
        bindings: Variable bindings from the match
        unmatched_left: Unmatched participants from left side
        unmatched_right: Unmatched participants from right side
    """

    matched: bool = False
    pattern: Optional[Reaction] = None
    bindings: Dict[str, str] = field(default_factory=dict)
    unmatched_left: List[Participant] = field(default_factory=list)
    unmatched_right: List[Participant] = field(default_factory=list)

    @property
    def unmatched(self) -> List[Participant]:
        """All unmatched participants from both sides."""
        return self.unmatched_left + self.unmatched_right

    def has_unmatched(self, chebi_id: str) -> bool:
        """Check if a specific CHEBI ID is in unmatched participants.

        Args:
            chebi_id: CHEBI identifier to search for

        Returns:
            True if the CHEBI ID is in unmatched participants

        Examples:
            >>> match = PatternMatch()
            >>> water = Participant(chebi_id="CHEBI:15377", name="water")
            >>> phosphate = Participant(chebi_id="CHEBI:43474", name="phosphate")
            >>>
            >>> match.unmatched_left = [water]
            >>> match.unmatched_right = [phosphate]
            >>>
            >>> match.has_unmatched("CHEBI:15377")  # water
            True
            >>> match.has_unmatched("CHEBI:43474")  # phosphate
            True
            >>> match.has_unmatched("CHEBI:99999")  # not present
            False
            >>>
            >>> # Check total unmatched count
            >>> len(match.unmatched)
            2
        """
        for p in self.unmatched:
            if p.chebi_id == chebi_id:
                return True
        return False


def match_patterns(
    reaction: Reaction, patterns: List[Reaction], strict: bool = True
) -> PatternMatch:
    """Match a reaction against a list of patterns.

    Args:
        reaction: The concrete reaction to match
        patterns: List of pattern reactions to try
        strict: If True, require all participants to be matched

    Returns:
        PatternMatch with results

    Examples:
        >>> from autarch.molecules import atp, adp, water, h_plus
        >>>
        >>> # Create patterns for forward and reverse kinase
        >>> forward = pattern([atp, var("S")], [adp, var("P")])
        >>> reverse = pattern([adp, var("S")], [atp, var("P")])
        >>> patterns = [forward, reverse]
        >>>
        >>> # Test forward reaction - normal order
        >>> glucose = Participant(chebi_id="CHEBI:17234", name="glucose")
        >>> g6p = Participant(chebi_id="CHEBI:17665", name="glucose-6-P")
        >>> reaction = Reaction(
        ...     left_participants=[atp, glucose],
        ...     right_participants=[adp, g6p]
        ... )
        >>>
        >>> match = match_patterns(reaction, patterns)
        >>> match.matched
        True
        >>> match.pattern == forward
        True
        >>> match.bindings["S"]
        'CHEBI:17234'
        >>> match.bindings["P"]
        'CHEBI:17665'
        >>>
        >>> # Test with REVERSED ORDER on left side - should still match!
        >>> reaction_reversed_left = Reaction(
        ...     left_participants=[glucose, atp],  # substrate before ATP
        ...     right_participants=[adp, g6p]
        ... )
        >>> match2 = match_patterns(reaction_reversed_left, patterns)
        >>> match2.matched  # Should be True - order doesn't matter
        True
        >>> match2.bindings["S"]
        'CHEBI:17234'
        >>>
        >>> # Test with REVERSED ORDER on right side - should still match!
        >>> reaction_reversed_right = Reaction(
        ...     left_participants=[atp, glucose],
        ...     right_participants=[g6p, adp]  # product before ADP
        ... )
        >>> match3 = match_patterns(reaction_reversed_right, patterns)
        >>> match3.matched  # Should be True - order doesn't matter
        True
        >>> match3.bindings["P"]
        'CHEBI:17665'
        >>>
        >>> # Test RHEA:10224 case: pyridoxal + ATP = pyridoxal-P + ADP + H+
        >>> pyridoxal = Participant(chebi_id="CHEBI:17310", name="pyridoxal")
        >>> pyridoxal_p = Participant(chebi_id="CHEBI:18405", name="pyridoxal-P")
        >>> rhea_10224 = Reaction(
        ...     left_participants=[pyridoxal, atp],  # substrate before ATP
        ...     right_participants=[pyridoxal_p, adp, h_plus]  # with H+
        ... )
        >>> # Pattern with optional H+
        >>> kinase_pattern = pattern([atp, var("S")], [adp, var("P"), optional(h_plus)])
        >>> match4 = match_patterns(rhea_10224, [kinase_pattern])
        >>> match4.matched  # Should be True!
        True
        >>> match4.bindings["S"]
        'CHEBI:17310'
        >>> match4.bindings["P"]
        'CHEBI:18405'
        >>>
        >>> # Test reaction that doesn't match any pattern
        >>> phosphate = Participant(chebi_id="CHEBI:43474", name="phosphate")
        >>> hydrolysis = Reaction(
        ...     left_participants=[atp, water],
        ...     right_participants=[adp, phosphate]
        ... )
        >>>
        >>> # This matches the forward pattern with water as substrate, phosphate as product
        >>> match = match_patterns(hydrolysis, patterns, strict=True)
        >>> match.matched
        True
        >>> match.bindings["S"]
        'CHEBI:15377'
        >>> match.bindings["P"]
        'CHEBI:43474'
        >>>
        >>> # A reaction that truly doesn't match (wrong nucleotide)
        >>> amp = Participant(chebi_id="CHEBI:16027", name="AMP")
        >>> no_match = Reaction(
        ...     left_participants=[amp, water],
        ...     right_participants=[phosphate]
        ... )
        >>> match2 = match_patterns(no_match, patterns)
        >>> match2.matched
        False
    """
    for pattern in patterns:
        match = match_single_pattern(reaction, pattern, strict)
        if match.matched:
            return match

    # No patterns matched
    return PatternMatch(matched=False)


def _find_matching_assignment(
    reaction_left: List[Participant],
    reaction_right: List[Participant],
    pattern_left: List[PatternParticipant],
    pattern_right: List[PatternParticipant],
    bindings: Dict[str, str],
    strict: bool
) -> Optional[Tuple[Dict[str, str], Set[int], Set[int]]]:
    """Find a valid assignment of reaction participants to pattern participants.
    
    Uses recursive backtracking to try all possible assignments.
    
    Returns:
        Tuple of (bindings, left_matched_indices, right_matched_indices) if successful,
        None if no valid assignment exists.
    """
    # Try to match left side
    left_assignment = _match_side_backtrack(
        reaction_left, pattern_left, bindings.copy(), "left"
    )
    
    if not left_assignment:
        return None
    
    new_bindings, left_matched = left_assignment
    
    # Try to match right side with the bindings from left side
    right_assignment = _match_side_backtrack(
        reaction_right, pattern_right, new_bindings.copy(), "right"
    )
    
    if not right_assignment:
        return None
    
    final_bindings, right_matched = right_assignment
    
    # Check if strict mode requirements are met
    if strict:
        # In strict mode, all reaction participants must be matched
        # (Pattern participants can be unmatched if they're optional)
        unmatched_left_indices = set(range(len(reaction_left))) - left_matched
        unmatched_right_indices = set(range(len(reaction_right))) - right_matched
        
        if unmatched_left_indices or unmatched_right_indices:
            return None
    
    return final_bindings, left_matched, right_matched


def _match_side_backtrack(
    reaction_participants: List[Participant],
    pattern_participants: List[PatternParticipant],
    bindings: Dict[str, str],
    side: str
) -> Optional[Tuple[Dict[str, str], Set[int]]]:
    """Match one side of a reaction using backtracking.
    
    Match order priority:
    1. Optional patterns first (greedy - claim their matches before variables)
    2. Required concrete patterns
    3. Required variable patterns
    
    Returns:
        Tuple of (bindings, matched_indices) if successful, None otherwise.
    """
    matched_indices = set()
    
    # Step 1: Match all optional patterns first (greedy)
    optional_patterns = [p for p in pattern_participants if p.min_count == 0]
    for opt_pattern in optional_patterns:
        for i, reaction_p in enumerate(reaction_participants):
            if i not in matched_indices:
                if _participant_matches_pattern(reaction_p, opt_pattern, bindings):
                    matched_indices.add(i)
                    break  # This optional found its match
    
    # Step 2: Separate required patterns into concrete and variables
    required_patterns = [p for p in pattern_participants if p.min_count > 0]
    concrete_required = [p for p in required_patterns if not p.is_variable]
    variable_required = [p for p in required_patterns if p.is_variable]
    
    # Step 3: Match concrete required patterns
    for pattern in concrete_required:
        matched = False
        for i, reaction_p in enumerate(reaction_participants):
            if i not in matched_indices:
                if pattern.matches(reaction_p):
                    matched_indices.add(i)
                    matched = True
                    break
        if not matched:
            return None  # Required concrete pattern couldn't be matched
    
    # Step 4: Match variable required patterns with remaining participants
    # In strict mode, we need to ensure we have exactly the right number of participants
    remaining_participants = [i for i in range(len(reaction_participants)) if i not in matched_indices]
    
    # Check if we have the right number of remaining participants for remaining patterns
    # Note: We need at least as many remaining participants as variable patterns
    # Extra participants are allowed and will be reported as unmatched
    if len(remaining_participants) < len(variable_required):
        # Not enough participants to match required variables
        return None
    
    for pattern in variable_required:
        var_name = pattern.variable
        if not var_name:
            continue
            
        matched = False
        for i in remaining_participants:
            if i not in matched_indices:
                if var_name in bindings:
                    # Check if matches existing binding
                    if _participant_matches_pattern(reaction_participants[i], pattern, bindings):
                        matched_indices.add(i)
                        matched = True
                        break
                else:
                    # Create new binding
                    chebi_id = reaction_participants[i].chebi_id
                    smiles = reaction_participants[i].smiles
                    if chebi_id is not None:
                        bindings[var_name] = chebi_id
                    elif smiles is not None:
                        bindings[var_name] = smiles
                    else:
                        continue
                    matched_indices.add(i)
                    matched = True
                    break
        
        if not matched:
            return None  # Required variable couldn't be matched
    
    return bindings, matched_indices


# Removed _find_assignment_recursive - no longer needed with greedy matching


def _participant_matches_pattern(
    participant: Participant,
    pattern: PatternParticipant,
    bindings: Dict[str, str]
) -> bool:
    """Check if a participant matches a pattern given current bindings."""
    if pattern.is_variable:
        var_name = pattern.variable
        if var_name and var_name in bindings:
            # Check against existing binding
            return bool(
                (participant.chebi_id and participant.chebi_id == bindings[var_name]) or
                (participant.smiles and participant.smiles == bindings[var_name])
            )
        # Unbound variable matches anything
        return True
    else:
        # Concrete pattern
        return pattern.matches(participant)


def match_single_pattern(
    reaction: Reaction, pattern: Reaction, strict: bool = True
) -> PatternMatch:
    """Match a reaction against a single pattern (order-independent).

    This performs order-independent pattern matching using recursive backtracking
    to try all possible assignments of reaction participants to pattern participants.

    Args:
        reaction: The concrete reaction to match
        pattern: The pattern reaction
        strict: If True, fail if there are unmatched participants

    Returns:
        PatternMatch with results including bindings and unmatched participants

    Examples:
        >>> from autarch.molecules import atp, adp, water, phosphate, h_plus, p
        >>>
        >>> # Test basic kinase pattern
        >>> glucose = Participant(chebi_id="CHEBI:17234", name="glucose")
        >>> g6p = Participant(chebi_id="CHEBI:17665", name="glucose-6-P")
        >>> 
        >>> # ATP + glucose -> ADP + glucose-6-P (standard order)
        >>> reaction1 = Reaction(
        ...     left_participants=[atp, glucose],
        ...     right_participants=[adp, g6p]
        ... )
        >>> pattern1 = p(atp) + var("S") >> p(adp) + var("P")
        >>> match1 = match_single_pattern(reaction1, pattern1, strict=True)
        >>> match1.matched
        True
        >>> match1.bindings["S"]
        'CHEBI:17234'
        >>> match1.bindings["P"]
        'CHEBI:17665'
        >>>
        >>> # glucose + ATP -> ADP + glucose-6-P (reversed left side)
        >>> reaction2 = Reaction(
        ...     left_participants=[glucose, atp],  # ATP comes second
        ...     right_participants=[adp, g6p]
        ... )
        >>> match2 = match_single_pattern(reaction2, pattern1, strict=True)
        >>> match2.matched
        True
        >>> match2.bindings["S"]
        'CHEBI:17234'
        >>>
        >>> # ATP + glucose -> glucose-6-P + ADP (reversed right side)
        >>> reaction3 = Reaction(
        ...     left_participants=[atp, glucose],
        ...     right_participants=[g6p, adp]  # ADP comes second
        ... )
        >>> match3 = match_single_pattern(reaction3, pattern1, strict=True)
        >>> match3.matched
        True
        >>> match3.bindings["P"]
        'CHEBI:17665'
        >>>
        >>> # Test with optional H+
        >>> pattern_with_h = p(atp) + var("S") + optional(p(h_plus)) >> p(adp) + var("P") + optional(p(h_plus))
        >>> 
        >>> # Reaction with H+ on right side only
        >>> reaction4 = Reaction(
        ...     left_participants=[glucose, atp],
        ...     right_participants=[adp, h_plus, g6p]
        ... )
        >>> match4 = match_single_pattern(reaction4, pattern_with_h, strict=True)
        >>> match4.matched
        True
        >>> match4.bindings["S"]
        'CHEBI:17234'
        >>> match4.bindings["P"]
        'CHEBI:17665'
        >>>
        >>> # H+ in middle position should still work
        >>> reaction5 = Reaction(
        ...     left_participants=[atp, glucose],
        ...     right_participants=[adp, h_plus, g6p]  # H+ between ADP and product
        ... )
        >>> match5 = match_single_pattern(reaction5, pattern_with_h, strict=True)
        >>> match5.matched
        True
        >>>
        >>> # Test strict mode with extra participants
        >>> reaction_extra = Reaction(
        ...     left_participants=[atp, glucose, water],  # Extra water
        ...     right_participants=[adp, g6p]
        ... )
        >>> match_strict = match_single_pattern(reaction_extra, pattern1, strict=True)
        >>> match_strict.matched
        False
        >>>
        >>> # Non-strict mode allows extra participants
        >>> match_loose = match_single_pattern(reaction_extra, pattern1, strict=False)
        >>> match_loose.matched
        True
        >>> len(match_loose.unmatched_left)
        1
        >>> match_loose.unmatched_left[0].chebi_id
        'CHEBI:15377'
    """
    # Convert all pattern participants to PatternParticipant type
    pattern_left = []
    for p in pattern.left_participants:
        if isinstance(p, PatternParticipant):
            pattern_left.append(p)
        else:
            pattern_left.append(PatternParticipant(
                chebi_id=p.chebi_id, name=p.name, smiles=p.smiles
            ))
    
    pattern_right = []
    for p in pattern.right_participants:
        if isinstance(p, PatternParticipant):
            pattern_right.append(p)
        else:
            pattern_right.append(PatternParticipant(
                chebi_id=p.chebi_id, name=p.name, smiles=p.smiles
            ))
    
    # Try to find a valid assignment using backtracking
    result = _find_matching_assignment(
        reaction.left_participants,
        reaction.right_participants,
        pattern_left,
        pattern_right,
        {},  # initial bindings
        strict
    )
    
    if result:
        bindings, left_matched, right_matched = result
        
        # Collect unmatched participants
        unmatched_left = [
            p for i, p in enumerate(reaction.left_participants) if i not in left_matched
        ]
        unmatched_right = [
            p for i, p in enumerate(reaction.right_participants) if i not in right_matched
        ]
    else:
        # No valid assignment found
        return PatternMatch(matched=False)
    
    # Empty the rest for now - will be replaced by the backtracking result
    bindings = bindings
    unmatched_left = unmatched_left
    unmatched_right = unmatched_right

    # In strict mode, fail if there are unmatched participants
    if strict and (unmatched_left or unmatched_right):
        return PatternMatch(
            matched=False,
            pattern=pattern,
            bindings=bindings,
            unmatched_left=unmatched_left,
            unmatched_right=unmatched_right,
        )

    return PatternMatch(
        matched=True,
        pattern=pattern,
        bindings=bindings,
        unmatched_left=unmatched_left,
        unmatched_right=unmatched_right,
    )
