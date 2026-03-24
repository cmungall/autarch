# Pattern Matching DSL for Reaction Classification

## Overview

This document describes a declarative pattern-matching approach for reaction classification, replacing imperative code with declarative patterns.

## Current Approach (Imperative)

```python
class Kinase(Transferase):
    def check_membership_impl(self, reaction: Reaction):
        # 150+ lines of conditional logic
        if has_atp_reactant and has_adp_product:
            if n_reactants != 2:
                return False
            # ... many more checks
```

## New Approach (Declarative)

```python
@pattern_classifier([
    ReactionPattern(
        name="Forward kinase",
        left=[
            ParticipantPattern(chebi_id="CHEBI:30616"),  # ATP
            ParticipantPattern(variable="X")             # substrate
        ],
        right=[
            ParticipantPattern(chebi_id="CHEBI:456216"),  # ADP
            ParticipantPattern(variable="X", modifications=["phosphorylated"]),
            ParticipantPattern(chebi_id="CHEBI:15378", min_count=0, max_count=1)  # H+
        ],
        require_no_water=True,
        forbidden_participants=["CHEBI:43474"]  # no free phosphate
    )
])
class Kinase(Transferase):
    pass  # Pattern does all the work!
```

## Key Design Elements

### 1. ParticipantPattern

Represents a pattern for matching reaction participants:

- **Exact match**: `chebi_id="CHEBI:30616"` (ATP)
- **Variable**: `variable="X"` (any molecule, can be cross-referenced)
- **Class match**: `chebi_class="CHEBI:36080"` (match any protein)
- **Optional**: `min_count=0, max_count=1` (e.g., H+ may or may not appear)
- **Exclusions**: `exclude_chebi_ids=["CHEBI:15377"]` (not water)
- **Modifications**: `modifications=["phosphorylated"]` (for tracking changes)

### 2. ReactionPattern

Defines a complete reaction pattern:

- **Left/right participants**: Lists of ParticipantPatterns
- **Global constraints**: `require_no_water`, `forbidden_participants`
- **Metadata**: Name and description for debugging

### 3. Variable Unification

Variables (like `X`) can appear on both sides and must unify to the same molecule:

```python
# X must be the same molecule on both sides
left=[ParticipantPattern(variable="X")]
right=[ParticipantPattern(variable="X")]
```

## Benefits

1. **Clarity**: Pattern is immediately understandable
2. **Maintainability**: Add new patterns without touching logic
3. **Testability**: Patterns can be tested in isolation
4. **Composability**: Patterns can be combined and reused
5. **Documentation**: Pattern IS the documentation

## Future Enhancements

### 1. Chemical Modifications

Currently, we track modifications as strings. Future work:

```python
@dataclass
class ChemicalModification:
    type: str  # "phosphorylation", "methylation", etc.
    group: str  # "PO3", "CH3", etc.
    mass_delta: float  # +79.966 for phosphorylation
    
    def verify(self, substrate_smiles: str, product_smiles: str) -> bool:
        """Check if product is modified substrate."""
        # Use RDKit to verify the modification
```

### 2. CHEBI Ontology Integration

```python
ParticipantPattern(
    chebi_class="CHEBI:36080",  # protein
    # Would match any protein: insulin, albumin, etc.
)
```

### 3. Complex Constraints

```python
@dataclass
class ReactionPattern:
    # ...
    validators: List[Callable] = field(default_factory=list)
    
    # Example: Check that product mass = substrate mass + 79.966 Da
    validators=[
        lambda reaction, bindings: check_phosphorylation(
            bindings["substrate"], 
            bindings["product"]
        )
    ]
```

### 4. Pattern Composition

```python
# Base pattern for any NTP kinase
NTP_KINASE = ReactionPattern(...)

# Specific patterns inherit/extend
ATP_KINASE = NTP_KINASE.with_participant(chebi_id="CHEBI:30616")
GTP_KINASE = NTP_KINASE.with_participant(chebi_id="CHEBI:37565")
```

### 5. Bidirectional Patterns

```python
@dataclass
class BidirectionalPattern(ReactionPattern):
    """Automatically matches forward and reverse."""
    
    def matches(self, reaction):
        # Try forward
        if super().matches(reaction):
            return True
        # Try reverse
        return self.reverse().matches(reaction)
```

### 6. SMARTS-like String DSL

For simple cases, a string DSL might be cleaner:

```python
KINASE_PATTERN = "CHEBI:30616 + $X >> CHEBI:456216 + $X~P + CHEBI:15378{0,1}"
```

Where:
- `$X` = variable
- `~P` = phosphorylated
- `{0,1}` = count range

### 7. Pattern Learning

Given a set of known kinase reactions, automatically learn patterns:

```python
def learn_pattern(reactions: List[Reaction], labels: List[bool]) -> ReactionPattern:
    """Learn a pattern from labeled examples."""
    # Find common participants, constraints, etc.
```

## Implementation Status

- ✅ Basic ParticipantPattern and ReactionPattern classes
- ✅ Simple variable unification
- ✅ Optional participants (min/max counts)
- ✅ Forbidden participants and water/phosphate constraints
- ✅ Pattern-based classifier decorator
- ✅ Test suite demonstrating functionality

## Next Steps

1. **Integrate with existing classifiers**: Gradually migrate from imperative to declarative
2. **CHEBI ontology support**: Enable class-based matching
3. **Chemical modification verification**: Use RDKit to verify transformations
4. **Performance optimization**: Pattern matching can be optimized with indexing
5. **Pattern validation**: Ensure patterns are chemically sensible

## Example Migration

### Before (150 lines):
```python
def check_membership_impl(self, reaction):
    # Complex conditional logic...
```

### After (10 lines):
```python
PATTERNS = [
    ReactionPattern(
        left=[ATP, substrate],
        right=[ADP, phospho_substrate, optional_H],
        constraints=[no_water, no_phosphate]
    )
]
```

## Conclusion

The pattern matching DSL provides a more maintainable, testable, and understandable way to define reaction classifiers. It separates the "what" (the pattern) from the "how" (the matching algorithm), making the codebase more modular and easier to extend.