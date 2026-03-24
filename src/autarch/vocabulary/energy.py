"""Energy-related molecules and currency metabolites.

This module contains molecules involved in energy metabolism,
including nucleotide phosphates and phosphate groups.
"""

from autarch.formula import molecule

# Adenosine phosphates
ATP = molecule("CHEBI:15422", name="ATP")
ADP = molecule("CHEBI:16761", name="ADP")
AMP = molecule("CHEBI:16027", name="AMP")

# Guanosine phosphates
GTP = molecule("CHEBI:15996", name="GTP")
GDP = molecule("CHEBI:17552", name="GDP")
GMP = molecule("CHEBI:17345", name="GMP")

# Phosphate groups
Pi = molecule("CHEBI:43474", name="Pi")  # Inorganic phosphate
PPi = molecule("CHEBI:33019", name="PPi")  # Pyrophosphate

# Additional energy molecules
CTP = molecule("CHEBI:17677", name="CTP")
CDP = molecule("CHEBI:17239", name="CDP")
CMP = molecule("CHEBI:17361", name="CMP")

UTP = molecule("CHEBI:15713", name="UTP")
UDP = molecule("CHEBI:17659", name="UDP")
UMP = molecule("CHEBI:16695", name="UMP")

# High-energy compounds
phosphoenolpyruvate = molecule("CHEBI:18021", name="PEP")
creatine_phosphate = molecule("CHEBI:17287", name="phosphocreatine")
acetyl_phosphate = molecule("CHEBI:22191", name="acetyl-P")

__all__ = [
    "ATP",
    "ADP",
    "AMP",
    "GTP",
    "GDP",
    "GMP",
    "CTP",
    "CDP",
    "CMP",
    "UTP",
    "UDP",
    "UMP",
    "Pi",
    "PPi",
    "phosphoenolpyruvate",
    "creatine_phosphate",
    "acetyl_phosphate",
]
