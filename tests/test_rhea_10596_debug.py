"""Debug RHEA:10596 specifically - this should clearly match the kinase pattern."""

from autarch.datamodel import Reaction, Participant
from autarch.ontology.kinase import Kinase
from autarch.pattern_dsl import match_patterns, match_single_pattern


def test_rhea_10596_should_match():
    """RHEA:10596: L-tyrosyl- + ATP = O-phospho-L-tyrosyl- + ADP + H(+)
    
    This is a textbook kinase reaction and should definitely match the pattern:
    ATP + substrate -> ADP + product + H+
    """
    
    # Create the reaction exactly as it appears in the evaluation
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),  # substrate
            Participant(chebi_id="CHEBI:30616", name="ATP"),                 # ATP
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),  # phosphorylated product
            Participant(chebi_id="CHEBI:456216", name="ADP"),                          # ADP
            Participant(chebi_id="CHEBI:15378", name="H(+)"),                          # proton
        ],
    )
    
    print("=== RHEA:10596 DETAILED DEBUG ===")
    print("Reaction: L-tyrosyl- + ATP = O-phospho-L-tyrosyl- + ADP + H(+)")
    print(f"Left participants: {[(p.chebi_id, p.name) for p in reaction.left_participants]}")
    print(f"Right participants: {[(p.chebi_id, p.name) for p in reaction.right_participants]}")
    
    kinase = Kinase()
    print("\nKinase patterns:")
    for i, pattern in enumerate(kinase.PATTERNS):
        print(f"  Pattern {i}:")
        left_details = []
        for p in pattern.left_participants:
            if p.variable:
                left_details.append(f"var('{p.variable}')")
            elif p.min_count == 0:
                left_details.append(f"optional({p.chebi_id or p.name})")
            else:
                left_details.append(f"{p.chebi_id or p.name}")
        
        right_details = []
        for p in pattern.right_participants:
            if p.variable:
                right_details.append(f"var('{p.variable}')")
            elif p.min_count == 0:
                right_details.append(f"optional({p.chebi_id or p.name})")
            else:
                right_details.append(f"{p.chebi_id or p.name}")
                
        print(f"    Left: {left_details}")
        print(f"    Right: {right_details}")
    
    # Test each pattern individually
    print("\nTesting individual patterns:")
    for i, pattern in enumerate(kinase.PATTERNS):
        match = match_single_pattern(reaction, pattern, strict=True)
        print(f"  Pattern {i}: matched={match.matched}")
        if match.matched:
            print(f"    Bindings: {match.bindings}")
        else:
            print(f"    Unmatched left: {[p.name for p in match.unmatched_left]}")
            print(f"    Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # Test overall pattern matching
    match = match_patterns(reaction, kinase.PATTERNS, strict=True)
    print("\nOverall pattern matching:")
    print(f"  Matched: {match.matched}")
    print(f"  Bindings: {match.bindings}")
    print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # Test the kinase classifier
    result = kinase.check_membership(reaction)
    print("\nKinase classification:")
    print(f"  is_member: {result.is_member}")
    print(f"  explanation: {result.explanation}")
    
    # This should absolutely match:
    # - Left: L-tyrosyl + ATP (matches: substrate + ATP)  
    # - Right: O-phospho-L-tyrosyl + ADP + H+ (matches: product + ADP + H+)
    # - Pattern: ATP + substrate + optional(H+) -> ADP + product + optional(H+)
    
    print("\n=== EXPECTED MATCHING ===")
    print("Left side should match:")
    print("  L-tyrosine residue -> variable 'substrate'")
    print("  ATP -> ATP pattern")
    print("Right side should match:")
    print("  L-tyrosine-phosphate residue -> variable 'product'") 
    print("  ADP -> ADP pattern")
    print("  H(+) -> optional H+ pattern")
    
    if not match.matched:
        print("\n!!! PATTERN MATCHING BUG DETECTED !!!")
        print("This reaction clearly fits the kinase pattern but failed to match.")
        
        # Let's debug the pattern matching step by step
        pattern_0 = kinase.PATTERNS[0]
        print("\nDebugging Pattern 0 matching:")
        print("Pattern left participants:")
        for j, p in enumerate(pattern_0.left_participants):
            print(f"  [{j}] {p.chebi_id} / {p.name} / var={p.variable} / min_count={p.min_count}")
        
        print("Pattern right participants:")
        for j, p in enumerate(pattern_0.right_participants):
            print(f"  [{j}] {p.chebi_id} / {p.name} / var={p.variable} / min_count={p.min_count}")
            
        print("Reaction left participants:")
        for j, p in enumerate(reaction.left_participants):
            print(f"  [{j}] {p.chebi_id} / {p.name}")
            
        print("Reaction right participants:")
        for j, p in enumerate(reaction.right_participants):
            print(f"  [{j}] {p.chebi_id} / {p.name}")
    
    assert match.matched, "RHEA:10596 should clearly match the kinase pattern - this is a bug!"


def test_simplified_rhea_10596():
    """Test RHEA:10596 with simplified participant ordering to isolate the issue."""
    
    # Try different orderings to see if that's the issue
    orderings = [
        # Original order
        {
            "name": "original",
            "left": [
                Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
                Participant(chebi_id="CHEBI:30616", name="ATP"),
            ],
            "right": [
                Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
            ],
        },
        # ATP first (like the pattern)
        {
            "name": "ATP_first",
            "left": [
                Participant(chebi_id="CHEBI:30616", name="ATP"),
                Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            ],
            "right": [
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
                Participant(chebi_id="CHEBI:15378", name="H(+)"),
            ],
        },
        # Without H+ to test optional matching
        {
            "name": "no_H_plus",
            "left": [
                Participant(chebi_id="CHEBI:30616", name="ATP"),
                Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            ],
            "right": [
                Participant(chebi_id="CHEBI:456216", name="ADP"),
                Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            ],
        },
    ]
    
    kinase = Kinase()
    
    for ordering in orderings:
        print(f"\n=== Testing {ordering['name']} ordering ===")
        
        reaction = Reaction(
            left_participants=ordering["left"],
            right_participants=ordering["right"],
        )
        
        match = match_patterns(reaction, kinase.PATTERNS, strict=True)
        print(f"Matched: {match.matched}")
        
        if match.matched:
            print(f"SUCCESS: {ordering['name']} ordering works!")
            print(f"Bindings: {match.bindings}")
        else:
            print(f"FAILED: {ordering['name']} ordering failed")
            print(f"Unmatched left: {[p.name for p in match.unmatched_left]}")
            print(f"Unmatched right: {[p.name for p in match.unmatched_right]}")


if __name__ == "__main__":
    test_rhea_10596_should_match()
    test_simplified_rhea_10596()