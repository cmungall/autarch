"""Common small molecules and gases.

This module contains simple molecules commonly found in biochemical reactions,
including water, gases, and small inorganic compounds.
"""

from autarch.formula import molecule

# Water and related
H2O = molecule("CHEBI:15377", name="H2O")
H = molecule("CHEBI:15378", name="H+")  # Proton
OH = molecule("CHEBI:16234", name="OH-")  # Hydroxide

# Common gases
O2 = molecule("CHEBI:15379", name="O2")
CO2 = molecule("CHEBI:16526", name="CO2")
CO = molecule("CHEBI:17245", name="CO")
N2 = molecule("CHEBI:17997", name="N2")
NO = molecule("CHEBI:16480", name="NO")
N2O = molecule("CHEBI:17045", name="N2O")
H2 = molecule("CHEBI:18276", name="H2")

# Nitrogen compounds
NH3 = molecule("CHEBI:16134", name="NH3")  # Ammonia
NH4 = molecule("CHEBI:28938", name="NH4+")  # Ammonium
NO2 = molecule("CHEBI:16301", name="NO2-")  # Nitrite
NO3 = molecule("CHEBI:17632", name="NO3-")  # Nitrate

# Sulfur compounds
H2S = molecule("CHEBI:29919", name="H2S")
SO2 = molecule("CHEBI:18422", name="SO2")
SO3 = molecule("CHEBI:29384", name="SO3")

# Reactive oxygen species
H2O2 = molecule("CHEBI:16240", name="H2O2")  # Hydrogen peroxide
O2_radical = molecule("CHEBI:18421", name="O2•-")  # Superoxide
OH_radical = molecule("CHEBI:29191", name="•OH")  # Hydroxyl radical

# Halides
HCl = molecule("CHEBI:17883", name="HCl")
HBr = molecule("CHEBI:47266", name="HBr")
HI = molecule("CHEBI:43451", name="HI")
HF = molecule("CHEBI:29228", name="HF")

__all__ = [
    "H2O",
    "H",
    "OH",
    "O2",
    "CO2",
    "CO",
    "N2",
    "NO",
    "N2O",
    "H2",
    "NH3",
    "NH4",
    "NO2",
    "NO3",
    "H2S",
    "SO2",
    "SO3",
    "H2O2",
    "O2_radical",
    "OH_radical",
    "HCl",
    "HBr",
    "HI",
    "HF",
]
