#!/usr/bin/env python
"""Demonstration of the chemical vocabulary package.

Shows how to use the organized vocabulary of common biochemical molecules.
"""

from autarch.formula import molecule  # For creating custom molecules

from autarch.vocabulary import (
    # Energy molecules
    ATP,
    ADP,
    AMP,
    GTP,
    GDP,
    Pi,
    PPi,
    # Cofactors
    NAD,
    NADH,
    NADP,
    NADPH,
    FAD,
    FADH2,
    CoA,
    # Common molecules
    H2O,
    H,
    O2,
    CO2,
    NH3,
    # Metabolites
    glucose,
    pyruvate,
    lactate,
    citrate,
    glutamate,
    Na,
    K,
    Ca,
    ALL_MOLECULES,
)

# Can also import from specific submodules
from autarch.vocabulary.energy import phosphoenolpyruvate
from autarch.vocabulary.cofactors import acetyl_CoA, SAM, SAH, biotin
from autarch.vocabulary.metabolites import (
    glucose_6P,
    fructose_6P,
    alpha_ketoglutarate,
    succinate,
    fumarate,
)


def demo_energy_reactions():
    """Demonstrate reactions with energy molecules."""
    print("\n=== ENERGY METABOLISM ===")

    # ATP hydrolysis
    hydrolysis = ATP + H2O >> ADP + Pi
    print(f"ATP hydrolysis: {hydrolysis}")

    # GTP hydrolysis
    gtp_hydrolysis = GTP + H2O >> GDP + Pi
    print(f"GTP hydrolysis: {gtp_hydrolysis}")

    # Pyrophosphate hydrolysis
    ppi_hydrolysis = PPi + H2O >> Pi * 2
    print(f"PPi hydrolysis: {ppi_hydrolysis}")

    # Substrate-level phosphorylation
    substrate_phosphorylation = phosphoenolpyruvate + ADP >> pyruvate + ATP
    print(f"PEP → pyruvate: {substrate_phosphorylation}")


def demo_redox_reactions():
    """Demonstrate redox reactions with cofactors."""
    print("\n=== REDOX REACTIONS ===")

    # NAD-dependent dehydrogenase
    ldh = lactate + NAD >> pyruvate + NADH + H
    print(f"Lactate dehydrogenase: {ldh}")

    # NADP-dependent reaction
    g6pd = (
        glucose_6P + NADP
        >> molecule("CHEBI:57925", name="6-phosphoglucono-lactone") + NADPH + H
    )
    print(f"G6P dehydrogenase: {g6pd}")

    # FAD-dependent reaction
    sdh = succinate + FAD >> fumarate + FADH2
    print(f"Succinate dehydrogenase: {sdh}")

    # Mixed cofactor system
    complex1 = NADH + H + ubiquinone >> NAD + ubiquinol
    print(f"Complex I (simplified): {complex1}")


def demo_metabolic_pathways():
    """Demonstrate central metabolic pathways."""
    print("\n=== METABOLIC PATHWAYS ===")

    # Glycolysis steps
    print("\nGlycolysis:")
    step1 = glucose + ATP >> glucose_6P + ADP + H
    print(f"  1. Hexokinase: {step1}")

    step2 = glucose_6P >> fructose_6P
    print(f"  2. Phosphoglucose isomerase: {step2}")

    # TCA cycle
    print("\nTCA Cycle:")
    tca1 = acetyl_CoA + oxaloacetate + H2O >> citrate + CoA + H
    print(f"  1. Citrate synthase: {tca1}")

    tca2 = alpha_ketoglutarate + NAD + CoA >> succinyl_CoA + NADH + CO2
    print(f"  2. α-KG dehydrogenase: {tca2}")

    # Amino acid metabolism
    print("\nAmino Acid Metabolism:")
    transamination = glutamate + oxaloacetate >> alpha_ketoglutarate + aspartate
    print(f"  Transamination: {transamination}")

    deamination = glutamate + NAD + H2O >> alpha_ketoglutarate + NADH + NH3
    print(f"  Glutamate dehydrogenase: {deamination}")


def demo_transport_reactions():
    """Demonstrate transport reactions with ions."""
    print("\n=== TRANSPORT REACTIONS ===")

    # Na/K pump
    na_k_pump = (Na["in"] * 3 + K["out"] * 2 + ATP) >> (
        Na["out"] * 3 + K["in"] * 2 + ADP + Pi
    )
    print(f"Na/K-ATPase: {na_k_pump}")

    # Ca pump
    ca_pump = (Ca["cytoplasm"] * 2 + ATP) >> (Ca["ER"] * 2 + ADP + Pi)
    print(f"Ca-ATPase: {ca_pump}")

    # Proton pump
    proton_pump = (H["cytoplasm"] + ATP) >> (H["extracellular"] + ADP + Pi)
    print(f"H+-ATPase: {proton_pump}")

    # Symporter
    na_glucose = (Na["out"] * 2 + glucose["out"]) >> (Na["in"] * 2 + glucose["in"])
    print(f"Na/glucose symporter: {na_glucose}")


def demo_cofactor_reactions():
    """Demonstrate reactions with various cofactors."""
    print("\n=== COFACTOR-DEPENDENT REACTIONS ===")

    # CoA reactions
    fa_activation = (
        palmitate + CoA + ATP
        >> molecule("CHEBI:15525", name="palmitoyl-CoA") + AMP + PPi
    )
    print(f"Fatty acid activation: {fa_activation}")

    # SAM methylation
    methylation = SAM + molecule("CHEBI:16040", name="substrate") >> SAH + molecule(
        "CHEBI:16041", name="methylated-substrate"
    )
    print(f"Methylation: {methylation}")

    # Biotin-dependent carboxylation
    carboxylation = acetyl_CoA + CO2 + ATP + biotin >> malonyl_CoA + ADP + Pi + biotin
    print(f"Acetyl-CoA carboxylase: {carboxylation}")

    # Metal-dependent reaction (simplified - Cu2 is imported at module level)
    from autarch.vocabulary.ions import Cu2
    from autarch.vocabulary.common import H2O2

    sod = O2_radical * 2 + H * 2 + Cu2 >> H2O2 + O2 + Cu2
    print(f"Superoxide dismutase: {sod}")


def demo_vocabulary_access():
    """Demonstrate accessing the vocabulary dictionary."""
    print("\n=== VOCABULARY ACCESS ===")

    print(f"\nTotal molecules in vocabulary: {len(ALL_MOLECULES)}")

    # Access molecule by name
    atp = ALL_MOLECULES["ATP"]
    print(f"ATP CHEBI ID: {atp.chebi_id}")
    print(f"ATP name: {atp.name}")

    # List categories
    print("\nMolecule categories:")
    print(
        f"  Energy molecules: {len([k for k in ALL_MOLECULES if k in ['ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'Pi', 'PPi']])}"
    )
    print(
        f"  Cofactors: {len([k for k in ALL_MOLECULES if 'NAD' in k or 'FAD' in k or 'CoA' in k])}"
    )
    print(
        f"  Amino acids: {len([k for k in ALL_MOLECULES if k in ['glutamate', 'glutamine', 'aspartate', 'asparagine']])}"
    )
    print(
        f"  Ions: {len([k for k in ALL_MOLECULES if k in ['Na', 'K', 'Ca', 'Mg', 'Fe2', 'Fe3']])}"
    )


def main():
    """Run all demonstrations."""
    print("=" * 60)
    print("  CHEMICAL VOCABULARY DEMONSTRATION")
    print("=" * 60)

    demo_energy_reactions()
    demo_redox_reactions()
    demo_metabolic_pathways()
    demo_transport_reactions()
    demo_cofactor_reactions()
    demo_vocabulary_access()

    print("\n" + "=" * 60)
    print("The vocabulary package provides:")
    print("  • 100+ pre-defined molecules with CHEBI IDs")
    print("  • Organized by category (energy, cofactors, metabolites, ions)")
    print("  • Easy import from submodules or main package")
    print("  • ALL_MOLECULES dictionary for programmatic access")
    print("=" * 60)


if __name__ == "__main__":
    # Import additional molecules for demos
    from autarch.vocabulary.metabolites import (
        fructose_6P,
        oxaloacetate,
        aspartate,
        palmitate,
    )
    from autarch.vocabulary.cofactors import (
        succinyl_CoA,
        malonyl_CoA,
        ubiquinone,
        ubiquinol,
    )
    from autarch.vocabulary.common import O2_radical

    main()
