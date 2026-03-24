"""Demo of the chemical reaction formula DSL."""

from autarch.formula import molecule, ATP, ADP, Pi, H2O, H, NAD, NADH, glucose, O2, CO2
from autarch.classifier import ReactionClassifier


def main():
    """Demonstrate the formula DSL for creating reactions."""

    print("=" * 60)
    print("Chemical Reaction Formula DSL Demo")
    print("=" * 60)

    # 1. Simple ATP hydrolysis
    print("\n1. ATP Hydrolysis (Hydrolase reaction):")
    atp_hydrolysis = ATP + H2O >> ADP + Pi
    print(f"   Formula: {atp_hydrolysis}")

    # Build and classify
    reaction = atp_hydrolysis.build()
    classifier = ReactionClassifier()
    results = classifier.classify(reaction)
    print("   Classification:")
    for name, result in results.items():
        if result.is_member:
            print(f"     ✓ {name}: {result.explanation}")

    # 2. Cellular respiration with stoichiometry
    print("\n2. Cellular Respiration (with stoichiometry):")
    respiration = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)
    print(f"   Formula: {respiration}")

    # 3. NAD reduction (Oxidoreductase)
    print("\n3. NAD Reduction (bidirectional):")
    nad_reduction = NAD + H | NADH
    print(f"   Formula: {nad_reduction}")

    # 4. ATP/ADP translocase (Transport reaction)
    print("\n4. ATP/ADP Translocase (with compartments):")
    translocase = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
    print(f"   Formula: {translocase}")

    reaction = translocase.build()
    print(f"   Is transport: {reaction.is_transport_reaction()}")
    transported = reaction.get_transported_molecules()
    for mol_id, from_loc, to_loc in transported:
        print(f"     {mol_id}: {from_loc} → {to_loc}")

    # 5. Kinase reaction
    print("\n5. Kinase Reaction (phosphorylation):")
    kinase = ATP + molecule("CHEBI:17234", name="glucose") >> ADP + molecule(
        "CHEBI:17665", name="glucose-6-P"
    )
    print(f"   Formula: {kinase}")

    # 6. Complex stoichiometry
    print("\n6. Complex Stoichiometry:")
    complex_rxn = (2 * H2O + O2 * 2) >> (4 * H) + (O2 * 2)
    print(f"   Formula: {complex_rxn}")

    # 7. Custom molecules
    print("\n7. Custom Molecules:")
    substrate = molecule("CHEBI:999", name="substrate")
    product = molecule("CHEBI:998", name="product")
    cofactor = molecule("CHEBI:997", name="cofactor")

    custom_rxn = (substrate + cofactor) >> (product + cofactor)
    print(f"   Formula: {custom_rxn}")

    # 8. Demonstrating operator precedence and grouping
    print("\n8. Operator Precedence:")
    # Multiplication binds tighter than addition
    rxn1 = ATP * 2 + H2O >> ADP * 2 + Pi * 2
    print(f"   Without parens: {rxn1}")
    print(f"   Left side has {len(rxn1.left)} participants")
    print(f"   ATP count: {rxn1.left.molecules[0].count}")

    # Note: You can't multiply a whole ParticipantList
    # (ATP + H2O) * 2 would fail because ParticipantList doesn't support multiplication
    # Instead, multiply individual molecules before combining:
    print("   Correct: ATP * 2 + H2O * 2 (multiply molecules individually)")

    print("\n" + "=" * 60)
    print("Key Features:")
    print("  • Use + to combine molecules")
    print("  • Use * for stoichiometry (2 * H2O or H2O * 2)")
    print("  • Use [] for compartments (ATP['cytoplasm'])")
    print("  • Use >> for left-to-right reactions")
    print("  • Use << for right-to-left reactions")
    print("  • Use | for bidirectional reactions")
    print("=" * 60)


if __name__ == "__main__":
    main()
