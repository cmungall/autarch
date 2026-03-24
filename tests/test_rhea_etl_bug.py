"""Test the RHEA ETL bug that causes stoichiometry duplication."""

from autarch.datamodel import Reaction, Participant
from autarch.ontology.kinase import Kinase


def test_rhea_10596_with_buggy_stoichiometry():
    """Test RHEA:10596 with the buggy 5x stoichiometry from the cache.
    
    This reproduces the exact bug causing the false negative.
    """
    
    # Recreate the EXACT data from the RHEA cache (with 5x duplication bug)
    reaction = Reaction(
        left_participants=[
            # 5 copies of L-tyrosine residue
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            # 5 copies of ATP
            Participant(chebi_id="CHEBI:30616", name="ATP"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
        ],
        right_participants=[
            # 5 copies of L-tyrosine-phosphate residue
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            # 5 copies of ADP
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            # 5 copies of H+
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
        ],
    )
    
    print("=== RHEA ETL BUG TEST ===")
    print(f"Left participants: {len(reaction.left_participants)} (should be 2, not 10)")
    print(f"Right participants: {len(reaction.right_participants)} (should be 3, not 15)")
    
    kinase = Kinase()
    
    print("\nKinase pattern expects:")
    pattern = kinase.PATTERNS[0]
    print(f"  Left: {len(pattern.left_participants)} participants")
    print(f"  Right: {len(pattern.right_participants)} participants")
    
    # Test classification
    result = kinase.check_membership(reaction)
    print("\nClassification result:")
    print(f"  is_member: {result.is_member}")
    print(f"  explanation: {result.explanation}")
    
    # This will fail because the pattern expects 3 left participants:
    # [ATP, substrate, optional(H+)] but we have 10 participants
    assert not result.is_member
    assert "No kinase pattern found" in result.explanation
    
    print("\n!!! CONFIRMED: RHEA ETL BUG CAUSES FALSE NEGATIVES !!!")
    print("The RHEA cache has duplicate participants causing stoichiometry mismatch.")


def test_correct_rhea_10596():
    """Test what RHEA:10596 should look like with correct stoichiometry."""
    
    # What the reaction SHOULD be (without duplication bug)
    reaction = Reaction(
        left_participants=[
            Participant(chebi_id="CHEBI:46858", name="L-tyrosine residue"),
            Participant(chebi_id="CHEBI:30616", name="ATP"),
        ],
        right_participants=[
            Participant(chebi_id="CHEBI:61978", name="L-tyrosine-phosphate residue"),
            Participant(chebi_id="CHEBI:456216", name="ADP"),
            Participant(chebi_id="CHEBI:15378", name="H(+)"),
        ],
    )
    
    print("\n=== CORRECT RHEA:10596 TEST ===")
    print(f"Left participants: {len(reaction.left_participants)}")
    print(f"Right participants: {len(reaction.right_participants)}")
    
    kinase = Kinase()
    result = kinase.check_membership(reaction)
    
    print("\nClassification result:")
    print(f"  is_member: {result.is_member}")
    print(f"  explanation: {result.explanation}")
    
    # This should work fine
    assert result.is_member
    assert "phosphoryl transfer" in result.explanation.lower()
    
    print("SUCCESS: Corrected RHEA:10596 is properly classified as kinase!")


if __name__ == "__main__":
    test_rhea_10596_with_buggy_stoichiometry()
    test_correct_rhea_10596()