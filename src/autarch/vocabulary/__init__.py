"""Chemical vocabulary for common molecules used in biochemical reactions.

This module provides pre-defined Molecule objects for commonly used chemicals
in biochemical reactions. Molecules are organized by category for easy access.

Examples:
    >>> from autarch.vocabulary import ATP, ADP, H2O, NAD, NADH
    >>> from autarch.vocabulary.energy import ATP, ADP, AMP
    >>> from autarch.vocabulary.cofactors import NAD, FAD
    >>> from autarch.vocabulary import ALL_MOLECULES
"""

from autarch.vocabulary.energy import (
    ATP,
    ADP,
    AMP,
    GTP,
    GDP,
    GMP,
    Pi,
    PPi,
)

from autarch.vocabulary.cofactors import (
    NAD,
    NADH,
    NADP,
    NADPH,
    FAD,
    FADH2,
    CoA,
    acetyl_CoA,
    THF,
    SAM,
    SAH,
)

from autarch.vocabulary.common import (
    H2O,
    H,
    OH,
    O2,
    CO2,
    NH3,
    H2O2,
    NO,
)

from autarch.vocabulary.metabolites import (
    glucose,
    fructose,
    galactose,
    pyruvate,
    lactate,
    citrate,
    succinate,
    fumarate,
    malate,
    glutamate,
    glutamine,
    aspartate,
    asparagine,
)

from autarch.vocabulary.ions import (
    Na,
    K,
    Ca,
    Mg,
    Cl,
    Fe2,
    Fe3,
    Cu2,
    Zn2,
    phosphate,
    sulfate,
    bicarbonate,
)

# Collect all molecules for easy access
ALL_MOLECULES = {
    # Energy
    "ATP": ATP,
    "ADP": ADP,
    "AMP": AMP,
    "GTP": GTP,
    "GDP": GDP,
    "GMP": GMP,
    "Pi": Pi,
    "PPi": PPi,
    # Cofactors
    "NAD": NAD,
    "NADH": NADH,
    "NADP": NADP,
    "NADPH": NADPH,
    "FAD": FAD,
    "FADH2": FADH2,
    "CoA": CoA,
    "acetyl_CoA": acetyl_CoA,
    "THF": THF,
    "SAM": SAM,
    "SAH": SAH,
    # Common
    "H2O": H2O,
    "H": H,
    "OH": OH,
    "O2": O2,
    "CO2": CO2,
    "NH3": NH3,
    "H2O2": H2O2,
    "NO": NO,
    # Metabolites
    "glucose": glucose,
    "fructose": fructose,
    "galactose": galactose,
    "pyruvate": pyruvate,
    "lactate": lactate,
    "citrate": citrate,
    "succinate": succinate,
    "fumarate": fumarate,
    "malate": malate,
    "glutamate": glutamate,
    "glutamine": glutamine,
    "aspartate": aspartate,
    "asparagine": asparagine,
    # Ions
    "Na": Na,
    "K": K,
    "Ca": Ca,
    "Mg": Mg,
    "Cl": Cl,
    "Fe2": Fe2,
    "Fe3": Fe3,
    "Cu2": Cu2,
    "Zn2": Zn2,
    "phosphate": phosphate,
    "sulfate": sulfate,
    "bicarbonate": bicarbonate,
}

__all__ = [
    # Energy molecules
    "ATP",
    "ADP",
    "AMP",
    "GTP",
    "GDP",
    "GMP",
    "Pi",
    "PPi",
    # Cofactors
    "NAD",
    "NADH",
    "NADP",
    "NADPH",
    "FAD",
    "FADH2",
    "CoA",
    "acetyl_CoA",
    "THF",
    "SAM",
    "SAH",
    # Common molecules
    "H2O",
    "H",
    "OH",
    "O2",
    "CO2",
    "NH3",
    "H2O2",
    "NO",
    # Metabolites
    "glucose",
    "fructose",
    "galactose",
    "pyruvate",
    "lactate",
    "citrate",
    "succinate",
    "fumarate",
    "malate",
    "glutamate",
    "glutamine",
    "aspartate",
    "asparagine",
    # Ions
    "Na",
    "K",
    "Ca",
    "Mg",
    "Cl",
    "Fe2",
    "Fe3",
    "Cu2",
    "Zn2",
    "phosphate",
    "sulfate",
    "bicarbonate",
    # Dictionary
    "ALL_MOLECULES",
]
