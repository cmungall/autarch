"""Debug stoichiometry issues in kinase pattern matching."""

from autarch.datamodel import Reaction, Participant
from autarch.ontology.kinase import Kinase
from autarch.pattern_dsl import match_patterns


def test_why_carbamoyl_phosphate_synthetase_fails():
    """Debug exactly why RHEA:10152 fails pattern matching."""
    
    # The failing reaction: hydrogencarbonate + NH4(+) + ATP = carbamoyl phosphate + ADP + H2O + H(+)
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:17544", name="hydrogencarbonate"),  # substrate 1
            Participant(chebi_id="CHEBI:28938", name="NH4(+)"),             # substrate 2
            Participant(chebi_id="CHEBI:30616", name="ATP"),                # ATP
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:58228", name="carbamoyl phosphate"), # product (incorporates both substrates)
            Participant(chebi_id="CHEBI:456216", name="ADP"),               # ADP
            Participant(chebi_id="CHEBI:15377", name="H2O"),                # water
            Participant(chebi_id="CHEBI:15378", name="H(+)"),               # proton
        ],
    )
    
    kinase = Kinase()
    
    print("=== CARBAMOYL PHOSPHATE SYNTHETASE DEBUG ===")
    print("Reaction participants:")
    print(f"  Left ({len(reaction.left_participants)}): {[p.name for p in reaction.left_participants]}")
    print(f"  Right ({len(reaction.right_participants)}): {[p.name for p in reaction.right_participants]}")
    
    print("\nKinase pattern expects:")
    for i, pattern in enumerate(kinase.PATTERNS):
        print(f"  Pattern {i} Left ({len(pattern.left_participants)}): {[p.name or f'?{p.variable}' for p in pattern.left_participants]}")
        print(f"  Pattern {i} Right ({len(pattern.right_participants)}): {[p.name or f'?{p.variable}' for p in pattern.right_participants]}")
    
    # Test pattern matching
    match = match_patterns(reaction, kinase.PATTERNS, strict=True)
    print("\nPattern matching (strict=True):")
    print(f"  Matched: {match.matched}")
    print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # Test with non-strict mode
    match_loose = match_patterns(reaction, kinase.PATTERNS, strict=False) 
    print("\nPattern matching (strict=False):")
    print(f"  Matched: {match_loose.matched}")
    print(f"  Unmatched left: {[p.name for p in match_loose.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match_loose.unmatched_right]}")
    
    # The issue: This is a legitimate kinase-like reaction but with complex stoichiometry
    # Pattern expects: ATP + substrate + H+ -> ADP + product + H+
    # Reality has: hydrogencarbonate + NH4+ + ATP -> carbamoyl_phosphate + ADP + H2O + H+
    #
    # In strict mode: Fails because there are extra participants (NH4+ left, H2O right)
    # In non-strict mode: Matches with unmatched participants reported

    assert not match.matched, "Expected to fail in strict mode due to stoichiometry mismatch"
    # Non-strict mode correctly matches and reports unmatched participants
    assert match_loose.matched, "Non-strict mode should match with unmatched participants"
    assert len(match_loose.unmatched_left) == 1  # NH4+ unmatched
    assert len(match_loose.unmatched_right) == 1  # H2O unmatched


def test_simple_kinase_that_works():
    """Test a simple kinase that works to contrast with the failing one."""
    
    # Simple kinase: D-glucosamine + ATP = D-glucosamine 6-phosphate + ADP + H(+) 
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:58723", name="D-glucosamine"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
            Participant(chebi_id="CHEBI:58725", name="D-glucosamine 6-phosphate"),
        ],
    )
    
    kinase = Kinase()
    
    print("\n=== SIMPLE GLUCOSAMINE KINASE DEBUG ===")
    print("Reaction participants:")
    print(f"  Left ({len(reaction.left_participants)}): {[p.name for p in reaction.left_participants]}")
    print(f"  Right ({len(reaction.right_participants)}): {[p.name for p in reaction.right_participants]}")
    
    match = match_patterns(reaction, kinase.PATTERNS, strict=True)
    print("\nPattern matching (strict=True):")
    print(f"  Matched: {match.matched}")
    print(f"  Bindings: {match.bindings}")
    print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # This works because:
    # - 2 left participants: D-glucosamine (substrate) + ATP 
    # - 3 right participants: ADP + H+ + D-glucosamine-6P (product)
    # - Pattern: ATP + substrate + optional(H+) -> ADP + product + optional(H+)
    # - Matching: ATP matches ATP, substrate matches D-glucosamine, product matches D-glucosamine-6P, H+ matches optional H+
    
    assert match.matched, "Simple kinase should match the pattern"


def test_atp_amp_kinase_missing_pattern():
    """Test ATP -> AMP kinase reactions that are not covered by current patterns."""
    
    # Example: ATP -> AMP + diphosphate (like pyruvate carboxylase)
    # RHEA:10756: pyruvate + phosphate + ATP = phosphoenolpyruvate + AMP + diphosphate + H(+)
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:15361", name="pyruvate"),
            Participant(chebi_id="CHEBI:43474", name="phosphate"), 
            Participant(chebi_id="CHEBI:30616", name="ATP"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:18021", name="phosphoenolpyruvate"),
            Participant(chebi_id="CHEBI:456215", name="AMP"),  # AMP not ADP!
            Participant(chebi_id="CHEBI:33019", name="diphosphate"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
        ],
    )
    
    kinase = Kinase()
    
    print("\n=== ATP->AMP KINASE DEBUG ===")
    print("Reaction participants:")
    print(f"  Left ({len(reaction.left_participants)}): {[p.name for p in reaction.left_participants]}")
    print(f"  Right ({len(reaction.right_participants)}): {[p.name for p in reaction.right_participants]}")
    
    match = match_patterns(reaction, kinase.PATTERNS, strict=True)
    print("\nPattern matching (strict=True):")
    print(f"  Matched: {match.matched}")
    print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # This fails because:
    # - Current patterns only cover ATP -> ADP
    # - This is ATP -> AMP + diphosphate  
    # - Need additional pattern for ATP -> AMP kinases
    
    assert not match.matched, "ATP->AMP kinases not covered by current patterns"


def test_ctp_kinase_missing_pattern():
    """Test CTP-based kinase reactions missing from current patterns."""
    
    # RHEA:10576: N-methylethanolamine phosphate + CTP + H(+) = CDP-N-methylethanolamine + diphosphate
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:57823", name="N-methylethanolamine phosphate"),
            Participant(chebi_id="CHEBI:37563", name="CTP"),  # CTP not ATP!
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:57828", name="CDP-N-methylethanolamine"),
            Participant(chebi_id="CHEBI:33019", name="diphosphate"),
        ],
    )
    
    kinase = Kinase()
    
    print("\n=== CTP KINASE DEBUG ===")
    print("Reaction participants:")
    print(f"  Left ({len(reaction.left_participants)}): {[p.name for p in reaction.left_participants]}")
    print(f"  Right ({len(reaction.right_participants)}): {[p.name for p in reaction.right_participants]}")
    
    match = match_patterns(reaction, kinase.PATTERNS, strict=True)
    print("\nPattern matching (strict=True):")
    print(f"  Matched: {match.matched}")
    print(f"  Unmatched left: {[p.name for p in match.unmatched_left]}")
    print(f"  Unmatched right: {[p.name for p in match.unmatched_right]}")
    
    # This fails because:
    # - Current patterns only cover ATP -> ADP and GTP -> GDP
    # - This is CTP -> CDP + diphosphate
    # - Need additional pattern for CTP -> CDP kinases
    
    assert not match.matched, "CTP->CDP kinases not covered by current patterns"