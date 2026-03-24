"""Ions and ionic compounds.

This module contains metal ions, anions, and other charged species
commonly found in biological systems.
"""

from autarch.formula import molecule

# Alkali metals (Group 1)
Li = molecule("CHEBI:49713", name="Li+")
Na = molecule("CHEBI:29101", name="Na+")
K = molecule("CHEBI:29103", name="K+")
Rb = molecule("CHEBI:33322", name="Rb+")
Cs = molecule("CHEBI:30440", name="Cs+")

# Alkaline earth metals (Group 2)
Mg = molecule("CHEBI:18420", name="Mg2+")
Ca = molecule("CHEBI:29108", name="Ca2+")
Sr = molecule("CHEBI:33324", name="Sr2+")
Ba = molecule("CHEBI:37134", name="Ba2+")

# Transition metals
Fe2 = molecule("CHEBI:29033", name="Fe2+")  # Ferrous
Fe3 = molecule("CHEBI:29034", name="Fe3+")  # Ferric
Cu = molecule("CHEBI:29036", name="Cu+")  # Cuprous
Cu2 = molecule("CHEBI:29035", name="Cu2+")  # Cupric
Zn2 = molecule("CHEBI:29105", name="Zn2+")
Mn2 = molecule("CHEBI:29035", name="Mn2+")
Co2 = molecule("CHEBI:48828", name="Co2+")
Ni2 = molecule("CHEBI:49786", name="Ni2+")
Mo = molecule("CHEBI:49868", name="Mo6+")  # Molybdate form

# Halides
F = molecule("CHEBI:17051", name="F-")
Cl = molecule("CHEBI:17996", name="Cl-")
Br = molecule("CHEBI:17203", name="Br-")
iodine = molecule("CHEBI:16382", name="I-")

# Phosphate species
phosphate = molecule("CHEBI:43474", name="PO4(3-)")  # Same as Pi
hydrogen_phosphate = molecule("CHEBI:43470", name="HPO4(2-)")
dihydrogen_phosphate = molecule("CHEBI:29888", name="H2PO4-")

# Sulfate species
sulfate = molecule("CHEBI:16189", name="SO4(2-)")
hydrogen_sulfate = molecule("CHEBI:29214", name="HSO4-")
sulfite = molecule("CHEBI:17359", name="SO3(2-)")
thiosulfate = molecule("CHEBI:16094", name="S2O3(2-)")

# Carbonate species
carbonate = molecule("CHEBI:17544", name="CO3(2-)")
bicarbonate = molecule("CHEBI:17544", name="HCO3-")

# Other anions
hydroxide = molecule("CHEBI:16234", name="OH-")
cyanide = molecule("CHEBI:17514", name="CN-")
acetate_ion = molecule("CHEBI:30089", name="CH3COO-")
oxalate = molecule("CHEBI:16995", name="C2O4(2-)")

# Complex ions
permanganate = molecule("CHEBI:25961", name="MnO4-")
dichromate = molecule("CHEBI:33141", name="Cr2O7(2-)")
ferricyanide = molecule("CHEBI:5020", name="Fe(CN)6(3-)")
ferrocyanide = molecule("CHEBI:5019", name="Fe(CN)6(4-)")

__all__ = [
    # Alkali metals
    "Li",
    "Na",
    "K",
    "Rb",
    "Cs",
    # Alkaline earth metals
    "Mg",
    "Ca",
    "Sr",
    "Ba",
    # Transition metals
    "Fe2",
    "Fe3",
    "Cu",
    "Cu2",
    "Zn2",
    "Mn2",
    "Co2",
    "Ni2",
    "Mo",
    # Halides
    "F",
    "Cl",
    "Br",
    "iodine",
    # Phosphate species
    "phosphate",
    "hydrogen_phosphate",
    "dihydrogen_phosphate",
    # Sulfate species
    "sulfate",
    "hydrogen_sulfate",
    "sulfite",
    "thiosulfate",
    # Carbonate species
    "carbonate",
    "bicarbonate",
    # Other anions
    "hydroxide",
    "cyanide",
    "acetate_ion",
    "oxalate",
    # Complex ions
    "permanganate",
    "dichromate",
    "ferricyanide",
    "ferrocyanide",
]
