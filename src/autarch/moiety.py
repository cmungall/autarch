"""SMARTS-based moiety detection with incremental caching.

Detects functional groups structurally instead of name-based matching.
Results are cached per (ChEBI ID, pattern_name) for efficiency.

Examples:
    >>> from autarch.moiety import Moiety, has_moiety, is_thioester
    >>> # Acetyl-CoA has a thioester bond
    >>> acetyl_coa_smiles = "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    >>> is_thioester(acetyl_coa_smiles)
    True
    >>> # Water has no thioester
    >>> is_thioester("O")
    False
    >>> # ATP has phosphate groups
    >>> atp_smiles = "Nc1ncnc2c1ncn2C3OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C3O"
    >>> has_moiety(atp_smiles, Moiety.PHOSPHATE)
    True
"""

import json
from enum import Enum
from pathlib import Path
from typing import Optional

from rdkit import Chem

CACHE_FILE = "cache/moiety_cache.json"
MOIETY_CACHE_VERSION = 2


class Moiety(Enum):
    """Common moieties with their SMARTS patterns.

    Each pattern is designed to detect a specific functional group
    in a molecule using RDKit SMARTS matching.
    """

    # Phosphate groups
    PHOSPHATE = "[P;X4](=[O;X1])([O;X2])([O;X2])[O;X2]"  # phosphate group
    PHOSPHOESTER = "[P;X4](=[O;X1])([O])([O])[O][C,N]"  # phosphate ester bond

    # Acyl/thioester groups
    THIOESTER = "[C;X3](=[O;X1])[S;X2]"  # C(=O)-S (CoA linkage)
    ACYL = "[C;X3](=[O;X1])[O,N,S;X2]"  # general acyl linkage (ester, amide, thioester)

    # Glycosidic bonds - sugar ring linked to non-ring
    O_GLYCOSIDE = "[C;R1][O;X2][C;!R1]"  # sugar ring O-linked to non-ring C
    N_GLYCOSIDE = "[C;R1][N;X3][C;!R1]"  # sugar ring N-linked to non-ring C

    # Other useful patterns
    CARBOXYL = "[C;X3](=[O;X1])[O;H1,X1-]"  # -COOH or -COO-
    HYDROXYL = "[O;H1][C;!$(C=O)]"  # -OH on non-carbonyl C
    AMINE = "[N;H2,H1;!$(N=*);!$(N#*)]"  # primary/secondary amine
    SULFATE = "[S;X4](=[O;X1])(=[O;X1])([O;X2])[O;X2]"  # sulfate group
    SULFAMATE = "[S;X4](=[O;X1])(=[O;X1])([N])[O]"  # S-N bond (sulfamate)

    # Aldehyde group (EC 1.2 substrates)
    ALDEHYDE = "[CX3H1,H2](=O)"  # R-CHO aldehyde, including formaldehyde
    KETONE = "[#6][CX3](=O)[#6]"  # R-CO-R ketone

    # Alpha-beta unsaturated carbonyl (enoyl, EC 1.3 substrates)
    # Note: [C;!R] excludes ring carbons to avoid matching NAD(P)H nicotinamide
    ENOYL = "[C;!R]=[C;!R][C](=O)[O,S,N]"  # C=C-C(=O)-X (enoyl-CoA/ACP pattern)
    ALPHA_BETA_UNSAT = "[C;!R]=[C;!R][C](=O)"  # General alpha-beta unsaturated carbonyl

    # Disulfide bond (not CH-OH oxidation)
    DISULFIDE = "[S;X2][S;X2]"  # S-S disulfide bond

    # Amino acid backbone (alpha carbon with amine and carboxyl)
    AMINO_ACID = "[NH2,NH3+][CH1]([*])[C](=O)[O]"  # alpha-amino acid pattern

    # ACP/protein linkage (phosphopantetheine)
    ACP_LINKAGE = "[C](=O)SCCNC(=O)CCNC"  # pantetheine arm pattern


# Global cache (loaded lazily)
_moiety_cache: Optional[dict[str, dict[str, bool]]] = None
_cache_modified: bool = False


def _load_cache() -> dict[str, dict[str, bool]]:
    """Load moiety cache from disk.

    Returns:
        Dict mapping ChEBI ID to dict of moiety name to bool.
    """
    global _moiety_cache
    if _moiety_cache is None:
        cache_path = Path(CACHE_FILE)
        if cache_path.exists():
            with open(cache_path) as f:
                cached_data = json.load(f)
            if (
                isinstance(cached_data, dict)
                and cached_data.get("__version__") == MOIETY_CACHE_VERSION
                and isinstance(cached_data.get("entries"), dict)
            ):
                _moiety_cache = cached_data["entries"]
            else:
                _moiety_cache = {}
        else:
            _moiety_cache = {}
    return _moiety_cache


def save_moiety_cache() -> None:
    """Save moiety cache to disk (call periodically or at end of processing)."""
    global _moiety_cache, _cache_modified
    if _moiety_cache and _cache_modified:
        cache_path = Path(CACHE_FILE)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(
                {"__version__": MOIETY_CACHE_VERSION, "entries": _moiety_cache},
                f,
                indent=2,
            )
        _cache_modified = False


def clear_moiety_cache() -> None:
    """Clear the in-memory cache (for testing)."""
    global _moiety_cache, _cache_modified
    _moiety_cache = None
    _cache_modified = False


def has_moiety(
    smiles: str,
    moiety: Moiety,
    chebi_id: Optional[str] = None
) -> bool:
    """Check if molecule contains a moiety using SMARTS.

    Args:
        smiles: SMILES string of the molecule
        moiety: Moiety enum value to check for
        chebi_id: Optional ChEBI ID for caching

    Returns:
        True if moiety is present

    Examples:
        >>> has_moiety("CC(=O)SC", Moiety.THIOESTER)  # thioester
        True
        >>> has_moiety("O", Moiety.THIOESTER)  # water
        False
        >>> has_moiety("OP(O)(O)=O", Moiety.PHOSPHATE)  # phosphoric acid
        True
    """
    global _moiety_cache, _cache_modified

    if not smiles:
        return False

    cache = _load_cache()

    # Check cache first (only if we have a ChEBI ID)
    if chebi_id and chebi_id in cache:
        if moiety.name in cache[chebi_id]:
            return cache[chebi_id][moiety.name]

    # Compute match
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    pattern = Chem.MolFromSmarts(moiety.value)
    if pattern is None:
        return False

    has_match = mol.HasSubstructMatch(pattern)

    # Cache result if we have a ChEBI ID
    if chebi_id:
        if chebi_id not in cache:
            cache[chebi_id] = {}
        cache[chebi_id][moiety.name] = has_match
        _cache_modified = True

    return has_match


def has_any_moiety(
    smiles: str,
    moieties: list[Moiety],
    chebi_id: Optional[str] = None
) -> bool:
    """Check if molecule contains any of the listed moieties.

    Args:
        smiles: SMILES string of the molecule
        moieties: List of Moiety enum values to check for
        chebi_id: Optional ChEBI ID for caching

    Returns:
        True if any moiety is present
    """
    return any(has_moiety(smiles, m, chebi_id) for m in moieties)


def has_all_moieties(
    smiles: str,
    moieties: list[Moiety],
    chebi_id: Optional[str] = None
) -> bool:
    """Check if molecule contains all of the listed moieties.

    Args:
        smiles: SMILES string of the molecule
        moieties: List of Moiety enum values to check for
        chebi_id: Optional ChEBI ID for caching

    Returns:
        True if all moieties are present
    """
    return all(has_moiety(smiles, m, chebi_id) for m in moieties)


# Convenience functions for common checks
def is_phosphorylated(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a phosphate group.

    Examples:
        >>> is_phosphorylated("OP(O)(O)=O")  # phosphoric acid
        True
        >>> is_phosphorylated("O")  # water
        False
    """
    return has_moiety(smiles, Moiety.PHOSPHATE, chebi_id)


def is_thioester(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a thioester bond (like CoA derivatives).

    Examples:
        >>> is_thioester("CC(=O)SC")  # methyl thioacetate
        True
        >>> is_thioester("CC(=O)OC")  # methyl acetate (ester, not thioester)
        False
    """
    return has_moiety(smiles, Moiety.THIOESTER, chebi_id)


def is_glycoside(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a glycosidic bond.

    Returns True for both O-glycosides and N-glycosides.
    """
    return has_any_moiety(smiles, [Moiety.O_GLYCOSIDE, Moiety.N_GLYCOSIDE], chebi_id)


def is_carboxylic_acid(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a carboxylic acid group."""
    return has_moiety(smiles, Moiety.CARBOXYL, chebi_id)


def is_sulfated(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a sulfate group."""
    return has_moiety(smiles, Moiety.SULFATE, chebi_id)


def is_sulfamate(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a sulfamate group (S-N bond).

    Sulfamates have an S-N bond where sulfur is bonded to nitrogen,
    e.g., R-SO2-NH2 or R-O-SO2-NH2.

    Examples:
        >>> is_sulfamate("NS(=O)(=O)OC1CCCCC1")  # cyclohexyl sulfamate
        True
        >>> is_sulfamate("O")  # water
        False
    """
    return has_moiety(smiles, Moiety.SULFAMATE, chebi_id)


def is_aldehyde(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains an aldehyde group (-CHO).

    Aldehydes are EC 1.2 substrates (aldehyde oxidoreductases).

    Examples:
        >>> is_aldehyde("CC=O")  # acetaldehyde
        True
        >>> is_aldehyde("C=O")  # formaldehyde
        True
        >>> is_aldehyde("CC(=O)C")  # acetone (ketone, not aldehyde)
        False
    """
    return has_moiety(smiles, Moiety.ALDEHYDE, chebi_id)


def is_ketone(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a ketone group (R-CO-R).

    Examples:
        >>> is_ketone("CC(=O)C")  # acetone
        True
        >>> is_ketone("CC=O")  # acetaldehyde
        False
    """
    return has_moiety(smiles, Moiety.KETONE, chebi_id)


def is_enoyl(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains an enoyl group (alpha-beta unsaturated acyl).

    Enoyl compounds are EC 1.3 substrates (CH-CH oxidoreductases).

    Examples:
        >>> is_enoyl("C=CC(=O)SC")  # simplified enoyl-thioester
        True
        >>> is_enoyl("CCC(=O)SC")  # saturated acyl-thioester
        False
    """
    return has_moiety(smiles, Moiety.ENOYL, chebi_id)


def has_disulfide(smiles: str, chebi_id: Optional[str] = None) -> bool:
    """Check if molecule contains a disulfide bond (S-S).

    Disulfide reactions are not CH-OH oxidation.
    """
    return has_moiety(smiles, Moiety.DISULFIDE, chebi_id)
