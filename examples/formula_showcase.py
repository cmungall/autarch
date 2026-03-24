#!/usr/bin/env python
"""Comprehensive showcase of the chemical reaction formula DSL.

This script demonstrates all features of the formula DSL including:
- Basic reaction creation
- Stoichiometry
- Compartments/locations
- Reaction directions
- Integration with the classifier system
"""

from autarch.formula import (
    molecule,
    ATP,
    ADP,
    Pi,
    H2O,
    H,
    NAD,
    NADH,
    FAD,
    FADH2,
    CoA,
    glucose,
    O2,
    CO2,
)
from autarch.classifier import ReactionClassifier


def section(title: str):
    """Print a section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def demo_basic_reactions():
    """Demonstrate basic reaction creation."""
    section("BASIC REACTIONS")

    # ATP hydrolysis
    print("\n1. ATP Hydrolysis:")
    hydrolysis = ATP + H2O >> ADP + Pi
    print(f"   Formula: {hydrolysis}")
    print(f"   Direction: {hydrolysis.direction}")

    # Build to Reaction object
    reaction = hydrolysis.build()
    print(f"   Left participants: {len(reaction.left_participants)}")
    print(f"   Right participants: {len(reaction.right_participants)}")

    # Classify
    classifier = ReactionClassifier()
    results = classifier.classify(reaction)
    print("   Classifications:")
    for name, result in results.items():
        if result.is_member:
            print(f"     ✓ {name}")


def demo_stoichiometry():
    """Demonstrate stoichiometry features."""
    section("STOICHIOMETRY")

    # Cellular respiration
    print("\n1. Cellular Respiration:")
    respiration = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)
    print(f"   Formula: {respiration}")

    reaction = respiration.build()
    print(f"   O2 count: {reaction.left_participants[1].count}")
    print(f"   CO2 count: {reaction.right_participants[0].count}")
    print(f"   H2O count: {reaction.right_participants[1].count}")

    # Water formation
    print("\n2. Water Formation:")
    # H2 doesn't exist as predefined, so create it
    H2 = molecule("CHEBI:18276", name="H2")
    water_formation = (H2 * 2) + O2 >> (H2O * 2)
    print(f"   Formula: {water_formation}")

    # Alternative syntax
    water_alt = (2 * H2) + O2 >> (2 * H2O)
    print(f"   Alternative: {water_alt}")

    # Mixed stoichiometry
    print("\n3. Complex Stoichiometry:")
    complex_rxn = (ATP * 3) + (H2O * 3) >> (ADP * 3) + (Pi * 3) + (H * 3)
    print(f"   Formula: {complex_rxn}")


def demo_compartments():
    """Demonstrate compartment/location features."""
    section("COMPARTMENTS AND TRANSPORT")

    # ATP/ADP translocase
    print("\n1. ATP/ADP Translocase:")
    translocase = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
    print(f"   Formula: {translocase}")

    reaction = translocase.build()
    print(f"   Is transport: {reaction.is_transport_reaction()}")

    transported = reaction.get_transported_molecules()
    print("   Transported molecules:")
    for mol_id, from_loc, to_loc in transported:
        mol_name = "ATP" if "15422" in mol_id else "ADP"
        print(f"     {mol_name}: {from_loc} → {to_loc}")

    # Proton pump
    print("\n2. Proton Pump:")
    proton_pump = (ATP["cytoplasm"] + H["cytoplasm"]) >> (
        ADP["cytoplasm"] + Pi["cytoplasm"] + H["extracellular"]
    )
    print(f"   Formula: {proton_pump}")

    # Glucose transporter
    print("\n3. Glucose Transporter:")
    glc_transport = glucose["extracellular"] >> glucose["cytoplasm"]
    print(f"   Formula: {glc_transport}")

    # Complex transport with stoichiometry
    print("\n4. Na+/glucose cotransporter:")
    Na = molecule("CHEBI:29101", name="Na+")
    cotransport = (Na["out"] * 2 + glucose["out"]) >> (Na["in"] * 2 + glucose["in"])
    print(f"   Formula: {cotransport}")


def demo_directions():
    """Demonstrate different reaction directions."""
    section("REACTION DIRECTIONS")

    # Left-to-right
    print("\n1. Left-to-right (>>):")
    forward = ATP + H2O >> ADP + Pi
    print(f"   Formula: {forward}")
    print(f"   Direction: {forward.direction}")

    # Right-to-left
    print("\n2. Right-to-left (<<):")
    reverse = ADP + Pi << ATP + H2O
    print(f"   Formula: {reverse}")
    print(f"   Direction: {reverse.direction}")
    print("   Note: This reverses the actual direction!")

    # Bidirectional
    print("\n3. Bidirectional (|):")
    equilibrium = NAD + H | NADH
    print(f"   Formula: {equilibrium}")
    print(f"   Direction: {equilibrium.direction}")

    # Another bidirectional example
    print("\n4. Lactate/Pyruvate equilibrium:")
    lactate = molecule("CHEBI:24996", name="lactate")
    pyruvate = molecule("CHEBI:15361", name="pyruvate")
    ldh_rxn = lactate + NAD | pyruvate + NADH + H
    print(f"   Formula: {ldh_rxn}")


def demo_custom_molecules():
    """Demonstrate creating custom molecules."""
    section("CUSTOM MOLECULES")

    # From CHEBI IDs
    print("\n1. From CHEBI IDs:")
    acetyl_coa = molecule("CHEBI:15351", name="acetyl-CoA")
    citrate = molecule("CHEBI:16947", name="citrate")
    oxaloacetate = molecule("CHEBI:16452", name="oxaloacetate")

    citrate_synthase = acetyl_coa + oxaloacetate + H2O >> citrate + CoA + H
    print(f"   Citrate synthase: {citrate_synthase}")

    # From SMILES
    print("\n2. From SMILES:")
    ethanol = molecule(smiles="CCO", name="ethanol")
    acetaldehyde = molecule(smiles="CC=O", name="acetaldehyde")

    alcohol_dehydrogenase = ethanol + NAD >> acetaldehyde + NADH + H
    print(f"   Alcohol dehydrogenase: {alcohol_dehydrogenase}")

    # Complex organic reaction
    print("\n3. Complex organic reaction:")
    benzene = molecule(smiles="c1ccccc1", name="benzene")
    nitrobenzene = molecule(smiles="c1ccc(cc1)[N+](=O)[O-]", name="nitrobenzene")
    HNO3 = molecule("CHEBI:48107", name="HNO3")

    nitration = benzene + HNO3 >> nitrobenzene + H2O
    print(f"   Nitration: {nitration}")


def demo_enzyme_classes():
    """Demonstrate how reactions map to enzyme classes."""
    section("ENZYME CLASSIFICATION")

    classifier = ReactionClassifier()

    # 1. Hydrolase
    print("\n1. Hydrolase Example:")
    hydrolase_rxn = ATP + H2O >> ADP + Pi
    print(f"   {hydrolase_rxn}")
    reaction = hydrolase_rxn.build()
    results = classifier.classify(reaction)
    print(f"   ✓ Hydrolase: {results['Hydrolase'].is_member}")

    # 2. Oxidoreductase
    print("\n2. Oxidoreductase Example:")
    oxidoreductase_rxn = NADH >> NAD + H
    print(f"   {oxidoreductase_rxn}")
    reaction = oxidoreductase_rxn.build()
    results = classifier.classify(reaction)
    print(f"   ✓ Oxidoreductase: {results['Oxidoreductase'].is_member}")

    # 3. Transferase (needs phosphate group transfer)
    print("\n3. Transferase Example:")
    substrate = molecule("CHEBI:999", name="substrate")
    substrate_P = molecule("CHEBI:998", name="substrate-P")
    transferase_rxn = ATP + substrate >> ADP + substrate_P
    print(f"   {transferase_rxn}")

    # 4. Lyase (C-C bond breaking)
    print("\n4. Lyase Example:")
    fructose_16_bp = molecule("CHEBI:16905", name="F-1,6-BP")
    g3p = molecule("CHEBI:17138", name="G3P")
    dhap = molecule("CHEBI:16108", name="DHAP")
    lyase_rxn = fructose_16_bp >> g3p + dhap
    print(f"   {lyase_rxn}")

    # 5. Isomerase
    print("\n5. Isomerase Example:")
    glucose_6p = molecule("CHEBI:14314", name="G6P")
    fructose_6p = molecule("CHEBI:15946", name="F6P")
    isomerase_rxn = glucose_6p >> fructose_6p
    print(f"   {isomerase_rxn}")

    # 6. Ligase
    print("\n6. Ligase Example:")
    amino_acid1 = molecule("CHEBI:33709", name="amino-acid-1")
    amino_acid2 = molecule("CHEBI:33708", name="amino-acid-2")
    dipeptide = molecule("CHEBI:46761", name="dipeptide")
    ligase_rxn = amino_acid1 + amino_acid2 + ATP >> dipeptide + ADP + Pi
    print(f"   {ligase_rxn}")


def demo_complex_pathways():
    """Demonstrate building metabolic pathways."""
    section("METABOLIC PATHWAYS")

    print("\n1. Glycolysis (first 3 steps):")

    # Step 1: Hexokinase
    glucose_6p = molecule("CHEBI:14314", name="G6P")
    step1 = glucose + ATP >> glucose_6p + ADP + H
    print(f"   Step 1 (Hexokinase): {step1}")

    # Step 2: Phosphoglucose isomerase
    fructose_6p = molecule("CHEBI:15946", name="F6P")
    step2 = glucose_6p >> fructose_6p
    print(f"   Step 2 (PGI): {step2}")

    # Step 3: Phosphofructokinase
    fructose_16bp = molecule("CHEBI:16905", name="F-1,6-BP")
    step3 = fructose_6p + ATP >> fructose_16bp + ADP + H
    print(f"   Step 3 (PFK): {step3}")

    print("\n2. TCA Cycle (selected reactions):")

    # Citrate synthase
    acetyl_coa = molecule("CHEBI:15351", name="acetyl-CoA")
    oxaloacetate = molecule("CHEBI:16452", name="OAA")
    citrate = molecule("CHEBI:16947", name="citrate")

    tca1 = acetyl_coa + oxaloacetate + H2O >> citrate + CoA + H
    print(f"   Citrate synthase: {tca1}")

    # Isocitrate dehydrogenase
    isocitrate = molecule("CHEBI:16219", name="isocitrate")
    alpha_kg = molecule("CHEBI:16810", name="α-KG")

    tca2 = isocitrate + NAD >> alpha_kg + NADH + CO2
    print(f"   Isocitrate DH: {tca2}")

    # Succinate dehydrogenase
    succinate = molecule("CHEBI:15741", name="succinate")
    fumarate = molecule("CHEBI:18012", name="fumarate")

    tca3 = succinate + FAD >> fumarate + FADH2
    print(f"   Succinate DH: {tca3}")


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("  CHEMICAL REACTION FORMULA DSL SHOWCASE")
    print("=" * 70)

    demo_basic_reactions()
    demo_stoichiometry()
    demo_compartments()
    demo_directions()
    demo_custom_molecules()
    demo_enzyme_classes()
    demo_complex_pathways()

    print("\n" + "=" * 70)
    print("  END OF SHOWCASE")
    print("=" * 70)
    print("\nThe formula DSL provides a natural way to express chemical reactions")
    print("using Python's operator overloading. Key operators:")
    print("  • + : combine molecules")
    print("  • * : stoichiometry")
    print("  • [] : compartments")
    print("  • >> : left-to-right")
    print("  • << : right-to-left")
    print("  • | : bidirectional")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
