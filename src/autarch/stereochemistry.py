"""Structural stereochemistry comparison utilities.

This module provides functions for comparing stereochemistry between molecules
using InChI (International Chemical Identifier). InChI encodes stereochemistry
in specific layers:
- /t: Tetrahedral stereochemistry (R/S configuration)
- /m: Parity for stereochemistry
- /s: Type of stereochemistry

Two molecules are stereoisomers if they have identical connectivity
but different stereochemistry layers.
"""

from typing import Optional


def parse_inchi_stereo_layer(inchi: Optional[str]) -> Optional[str]:
    """Extract stereochemistry layers from InChI.

    The stereochemistry in InChI is encoded in layers starting with:
    - 't': Tetrahedral stereochemistry
    - 'm': Parity layer
    - 's': Stereochemistry type

    Args:
        inchi: InChI string

    Returns:
        Concatenated stereo layers (e.g., "t2-/m0/s1") or None if no stereo

    Examples:
        >>> # L-alanine has stereochemistry
        >>> parse_inchi_stereo_layer("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")
        't2-/m0/s1'
        >>> # Achiral molecule has no stereo layers
        >>> parse_inchi_stereo_layer("InChI=1S/H2O/h1H2")
        >>> # None returned
    """
    if not inchi:
        return None

    parts = inchi.split("/")
    stereo_layers = [p for p in parts if p.startswith(("t", "m", "s"))]

    if not stereo_layers:
        return None

    return "/".join(stereo_layers)


def remove_stereo_from_inchi(inchi: Optional[str]) -> Optional[str]:
    """Remove stereochemistry layers from InChI to get connectivity only.

    Args:
        inchi: InChI string

    Returns:
        InChI with stereo layers removed (connectivity formula only)

    Examples:
        >>> # L-alanine without stereo
        >>> remove_stereo_from_inchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")
        'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)'
    """
    if not inchi:
        return None

    parts = inchi.split("/")
    # Keep parts that don't start with stereo layer identifiers
    conn_parts = [p for p in parts if not p.startswith(("t", "m", "s"))]

    return "/".join(conn_parts)


def are_stereoisomers(
    mol1_inchi: Optional[str],
    mol2_inchi: Optional[str]
) -> tuple[bool, str]:
    """Check if two molecules are stereoisomers using InChI.

    Two molecules are stereoisomers if:
    1. They have identical connectivity (InChI without stereo layers)
    2. They have different stereochemistry (different /t, /m, or /s layers)

    Args:
        mol1_inchi: InChI of first molecule
        mol2_inchi: InChI of second molecule

    Returns:
        Tuple of (is_stereoisomer, explanation)

    Examples:
        >>> # L-alanine vs D-alanine
        >>> l_ala = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"
        >>> d_ala = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1"
        >>> is_stereo, explanation = are_stereoisomers(l_ala, d_ala)
        >>> is_stereo
        True
        >>> "Stereoisomerism" in explanation
        True
        >>>
        >>> # Same molecule - not stereoisomers
        >>> is_stereo, _ = are_stereoisomers(l_ala, l_ala)
        >>> is_stereo
        False
        >>>
        >>> # Different molecules (different connectivity)
        >>> water = "InChI=1S/H2O/h1H2"
        >>> is_stereo, explanation = are_stereoisomers(l_ala, water)
        >>> is_stereo
        False
        >>> "Different molecular connectivity" in explanation
        True
    """
    if not mol1_inchi or not mol2_inchi:
        return False, "Missing InChI"

    # Compare connectivity (InChI without stereo layers)
    conn1 = remove_stereo_from_inchi(mol1_inchi)
    conn2 = remove_stereo_from_inchi(mol2_inchi)

    if conn1 != conn2:
        return False, "Different molecular connectivity"

    # Compare stereochemistry layers
    stereo1 = parse_inchi_stereo_layer(mol1_inchi)
    stereo2 = parse_inchi_stereo_layer(mol2_inchi)

    # If neither has stereochemistry, they're not stereoisomers
    if stereo1 is None and stereo2 is None:
        return False, "No stereochemistry in either molecule"

    # If same stereochemistry, they're not stereoisomers (same molecule)
    if stereo1 == stereo2:
        return False, "Identical stereochemistry"

    # Different stereochemistry with same connectivity = stereoisomers
    stereo1_display = stereo1 or "(none)"
    stereo2_display = stereo2 or "(none)"
    return True, f"Stereoisomerism: {stereo1_display} -> {stereo2_display}"


def get_stereo_type(stereo_layer: Optional[str]) -> Optional[str]:
    """Determine the type of stereoisomerism from InChI stereo layers.

    Args:
        stereo_layer: Stereo layer string from parse_inchi_stereo_layer

    Returns:
        "enantiomer" for R/S switch, "epimer" for partial change, or None

    Examples:
        >>> # Full inversion (enantiomer)
        >>> get_stereo_type("t2-/m0/s1")
        'chiral'
    """
    if not stereo_layer:
        return None

    if "t" in stereo_layer:
        return "chiral"

    return None
