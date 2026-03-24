"""ChEBI ID normalization for RHEA reactions.

This module provides functions to normalize and recover ChEBI IDs from:
1. EQUATION field in rhea-reactions.txt.gz (most reliable)
2. Name-based lookup for common molecules

Examples:
    >>> from autarch.etl.chebi_normalization import (
    ...     parse_equation_chebi_ids,
    ...     normalize_chebi_by_name,
    ...     COMMON_MOLECULE_CHEBI,
    ... )
    >>> parse_equation_chebi_ids("CHEBI:16459 + CHEBI:15377 = CHEBI:31011")
    (['CHEBI:16459', 'CHEBI:15377'], ['CHEBI:31011'])
    >>> normalize_chebi_by_name("H(+)")
    'CHEBI:15378'
"""

import re
from typing import Dict, List, Optional, Tuple


# Common molecules with their canonical ChEBI IDs
# These are molecules that frequently appear in RHEA but may not have
# ChEBI IDs assigned through SMILES lookup
COMMON_MOLECULE_CHEBI: Dict[str, str] = {
    # Protons and hydrons
    "h(+)": "CHEBI:15378",
    "h+": "CHEBI:15378",
    "proton": "CHEBI:15378",
    "hydron": "CHEBI:15378",
    "hydrogen ion": "CHEBI:15378",
    # Water
    "h2o": "CHEBI:15377",
    "water": "CHEBI:15377",
    "h(2)o": "CHEBI:15377",
    # Phosphate species
    "phosphate": "CHEBI:43474",
    "pi": "CHEBI:43474",
    "orthophosphate": "CHEBI:43474",
    "inorganic phosphate": "CHEBI:43474",
    "diphosphate": "CHEBI:33019",
    "ppi": "CHEBI:33019",
    "pyrophosphate": "CHEBI:33019",
    # Nucleotides - ATP family
    "atp": "CHEBI:30616",
    "adenosine triphosphate": "CHEBI:30616",
    "atp(4-)": "CHEBI:30616",
    "adp": "CHEBI:456216",
    "adenosine diphosphate": "CHEBI:456216",
    "adp(3-)": "CHEBI:456216",
    "amp": "CHEBI:456215",
    "adenosine monophosphate": "CHEBI:456215",
    "amp(2-)": "CHEBI:456215",
    # Nucleotides - GTP family
    "gtp": "CHEBI:37565",
    "guanosine triphosphate": "CHEBI:37565",
    "gtp(4-)": "CHEBI:37565",
    "gdp": "CHEBI:58189",
    "guanosine diphosphate": "CHEBI:58189",
    "gdp(3-)": "CHEBI:58189",
    "gmp": "CHEBI:58115",
    "guanosine monophosphate": "CHEBI:58115",
    # Nucleotides - UTP family
    "utp": "CHEBI:46398",
    "uridine triphosphate": "CHEBI:46398",
    "udp": "CHEBI:58223",
    "uridine diphosphate": "CHEBI:58223",
    "ump": "CHEBI:57865",
    "uridine monophosphate": "CHEBI:57865",
    # Nucleotides - CTP family
    "ctp": "CHEBI:37563",
    "cytidine triphosphate": "CHEBI:37563",
    "cdp": "CHEBI:58069",
    "cytidine diphosphate": "CHEBI:58069",
    "cmp": "CHEBI:60377",
    "cytidine monophosphate": "CHEBI:60377",
    # NAD/NADP cofactors
    "nad(+)": "CHEBI:57540",
    "nad+": "CHEBI:57540",
    "nadh": "CHEBI:57945",
    "nadp(+)": "CHEBI:58349",
    "nadp+": "CHEBI:58349",
    "nadph": "CHEBI:57783",
    # FAD/FMN cofactors
    "fad": "CHEBI:57692",
    "fadh2": "CHEBI:58307",
    "fmn": "CHEBI:58210",
    "fmnh2": "CHEBI:57618",
    # CoA and derivatives
    "coa": "CHEBI:57287",
    "coenzyme a": "CHEBI:57287",
    "coa-sh": "CHEBI:57287",
    "acetyl-coa": "CHEBI:57288",
    "malonyl-coa": "CHEBI:57384",
    # Gases
    "o2": "CHEBI:15379",
    "oxygen": "CHEBI:15379",
    "molecular oxygen": "CHEBI:15379",
    "dioxygen": "CHEBI:15379",
    "co2": "CHEBI:16526",
    "carbon dioxide": "CHEBI:16526",
    "h2": "CHEBI:18276",
    "hydrogen": "CHEBI:18276",
    "molecular hydrogen": "CHEBI:18276",
    "dihydrogen": "CHEBI:18276",
    "nh3": "CHEBI:16134",
    "ammonia": "CHEBI:16134",
    "nh4(+)": "CHEBI:28938",
    "nh4+": "CHEBI:28938",
    "ammonium": "CHEBI:28938",
    "ammonium ion": "CHEBI:28938",
    # Sulfur compounds
    "sulfur": "CHEBI:17909",
    "s": "CHEBI:17909",
    "sulfur atom": "CHEBI:17909",
    "hydrogen sulfide": "CHEBI:16136",
    "h2s": "CHEBI:16136",
    "sulfide": "CHEBI:16136",
    "sulfate": "CHEBI:16189",
    "so4(2-)": "CHEBI:16189",
    # Common metabolites
    "acetate": "CHEBI:30089",
    "pyruvate": "CHEBI:15361",
    "lactate": "CHEBI:24996",
    "succinate": "CHEBI:30031",
    "fumarate": "CHEBI:29806",
    "malate": "CHEBI:15589",
    "oxaloacetate": "CHEBI:16452",
    "citrate": "CHEBI:16947",
    "alpha-ketoglutarate": "CHEBI:16810",
    "2-oxoglutarate": "CHEBI:16810",
    # Amino acids (common forms)
    "l-glutamate": "CHEBI:29985",
    "glutamate": "CHEBI:29985",
    "l-glutamine": "CHEBI:58359",
    "glutamine": "CHEBI:58359",
    "l-aspartate": "CHEBI:29991",
    "aspartate": "CHEBI:29991",
    "l-asparagine": "CHEBI:58048",
    "asparagine": "CHEBI:58048",
    "glycine": "CHEBI:57305",
    # Sugars
    "glucose": "CHEBI:17234",
    "d-glucose": "CHEBI:17234",
    "fructose": "CHEBI:15903",
    "d-fructose": "CHEBI:15903",
    "galactose": "CHEBI:16551",
    "d-galactose": "CHEBI:16551",
    # Other common
    "formate": "CHEBI:15740",
    "formaldehyde": "CHEBI:16842",
    "methanol": "CHEBI:17790",
    "ethanol": "CHEBI:16236",
    "acetaldehyde": "CHEBI:15343",
    "hydrogen peroxide": "CHEBI:16240",
    "h2o2": "CHEBI:16240",
    "superoxide": "CHEBI:18421",
}


def normalize_chebi_by_name(name: str) -> Optional[str]:
    """Look up ChEBI ID for a molecule by its common name.

    Args:
        name: Molecule name (case-insensitive)

    Returns:
        ChEBI ID if found, None otherwise

    Examples:
        >>> normalize_chebi_by_name("H(+)")
        'CHEBI:15378'
        >>> normalize_chebi_by_name("water")
        'CHEBI:15377'
        >>> normalize_chebi_by_name("ATP")
        'CHEBI:30616'
        >>> normalize_chebi_by_name("unknown_molecule") is None
        True
    """
    if not name:
        return None
    normalized = name.lower().strip()
    return COMMON_MOLECULE_CHEBI.get(normalized)


def parse_equation_chebi_ids(equation: str) -> Tuple[List[str], List[str]]:
    """Parse ChEBI IDs from a RHEA EQUATION line.

    The EQUATION field contains ChEBI IDs in the format:
    - Left side + Right side separated by =, =>, or <=>
    - Participants separated by +
    - Stoichiometry prefix like "2 CHEBI:xxxxx"
    - Multiple instances as "CHEBI:xxx,CHEBI:xxx"

    Args:
        equation: EQUATION string from RHEA (without "EQUATION" prefix)

    Returns:
        Tuple of (left_chebi_ids, right_chebi_ids)

    Examples:
        >>> parse_equation_chebi_ids("CHEBI:16459 + CHEBI:15377 = CHEBI:31011")
        (['CHEBI:16459', 'CHEBI:15377'], ['CHEBI:31011'])
        >>> parse_equation_chebi_ids("2 CHEBI:15378 + CHEBI:30616 => CHEBI:456216")
        (['CHEBI:15378', 'CHEBI:15378', 'CHEBI:30616'], ['CHEBI:456216'])
        >>> parse_equation_chebi_ids("CHEBI:29950,CHEBI:29950 = CHEBI:50058")
        (['CHEBI:29950', 'CHEBI:29950'], ['CHEBI:50058'])
        >>> parse_equation_chebi_ids("")
        ([], [])
    """
    if not equation or not equation.strip():
        return [], []

    # Split by reaction arrow (handle <=>, =>, =)
    # Order matters: check <=> and => before =
    if "<=>" in equation:
        parts = equation.split("<=>")
    elif "=>" in equation:
        parts = equation.split("=>")
    elif "=" in equation:
        parts = equation.split("=")
    else:
        return [], []

    if len(parts) != 2:
        return [], []

    left_str, right_str = parts

    def extract_chebi_ids(side: str) -> List[str]:
        """Extract ChEBI IDs from one side of the equation."""
        chebi_ids = []
        # Split by +
        participants = side.split("+")
        for participant in participants:
            participant = participant.strip()
            if not participant:
                continue

            # Check for stoichiometry prefix (e.g., "2 CHEBI:15378" or "2n CHEBI:xxx")
            stoich_match = re.match(r"^(\d+)\s+(CHEBI:\d+)", participant)
            if stoich_match:
                count = int(stoich_match.group(1))
                chebi_id = stoich_match.group(2)
                chebi_ids.extend([chebi_id] * count)
                continue

            # Check for variable stoichiometry (e.g., "n CHEBI:xxx", "2n CHEBI:xxx")
            var_stoich_match = re.match(r"^(\d*n)\s+(CHEBI:\d+)", participant)
            if var_stoich_match:
                # For variable stoichiometry, just include once
                chebi_id = var_stoich_match.group(2)
                chebi_ids.append(chebi_id)
                continue

            # Check for comma-separated ChEBI IDs (e.g., "CHEBI:29950,CHEBI:29950")
            comma_ids = re.findall(r"CHEBI:\d+", participant)
            if comma_ids:
                chebi_ids.extend(comma_ids)

        return chebi_ids

    return extract_chebi_ids(left_str), extract_chebi_ids(right_str)


def parse_equation_to_position_map(
    equation: str,
) -> Tuple[Dict[int, str], Dict[int, str]]:
    """Parse EQUATION and return position-indexed ChEBI IDs.

    This is useful for matching ChEBI IDs to participants by position.

    Args:
        equation: EQUATION string from RHEA

    Returns:
        Tuple of (left_position_map, right_position_map)
        where each map is {position: chebi_id}

    Examples:
        >>> left, right = parse_equation_to_position_map(
        ...     "CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938"
        ... )
        >>> left
        {0: 'CHEBI:16459', 1: 'CHEBI:15377'}
        >>> right
        {0: 'CHEBI:31011', 1: 'CHEBI:28938'}
    """
    left_ids, right_ids = parse_equation_chebi_ids(equation)

    left_map = {i: chebi_id for i, chebi_id in enumerate(left_ids)}
    right_map = {i: chebi_id for i, chebi_id in enumerate(right_ids)}

    return left_map, right_map


class EquationChEBILookup:
    """Lookup ChEBI IDs from RHEA EQUATION data.

    This class loads EQUATION data from rhea-reactions.txt.gz and provides
    fast lookup of ChEBI IDs by RHEA ID.

    Examples:
        >>> lookup = EquationChEBILookup()
        >>> lookup.load_from_file("cache/rhea_tsv/rhea-reactions.txt.gz")  # doctest: +SKIP
        >>> left, right = lookup.get_chebi_ids("10000")  # doctest: +SKIP
    """

    def __init__(self):
        """Initialize empty lookup."""
        self._equations: Dict[str, str] = {}
        self._parsed_cache: Dict[str, Tuple[List[str], List[str]]] = {}

    def load_from_file(self, filepath: str) -> int:
        """Load EQUATION data from rhea-reactions.txt.gz.

        Args:
            filepath: Path to rhea-reactions.txt.gz

        Returns:
            Number of equations loaded

        Examples:
            >>> import tempfile
            >>> import gzip
            >>> lookup = EquationChEBILookup()
            >>> with tempfile.NamedTemporaryFile(suffix='.gz', delete=False) as f:
            ...     content = b'''ENTRY       RHEA:10000
            ... DEFINITION  test reaction
            ... EQUATION    CHEBI:12345 = CHEBI:67890
            ... ///
            ... '''
            ...     with gzip.open(f.name, 'wb') as gz:
            ...         _ = gz.write(content)
            ...     count = lookup.load_from_file(f.name)
            >>> count
            1
            >>> lookup.get_equation("10000")
            'CHEBI:12345 = CHEBI:67890'
        """
        import gzip

        current_rhea_id = None
        count = 0

        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip()

                if line.startswith("ENTRY"):
                    # Parse RHEA ID: "ENTRY       RHEA:10000"
                    match = re.match(r"ENTRY\s+RHEA:(\d+)", line)
                    if match:
                        current_rhea_id = match.group(1)

                elif line.startswith("EQUATION") and current_rhea_id:
                    # Parse equation: "EQUATION    CHEBI:xxx + CHEBI:yyy = ..."
                    equation = line[12:].strip()  # Skip "EQUATION    "
                    self._equations[current_rhea_id] = equation
                    count += 1

                elif line.startswith("///"):
                    current_rhea_id = None

        return count

    def get_equation(self, rhea_id: str) -> Optional[str]:
        """Get raw EQUATION string for a RHEA ID.

        Args:
            rhea_id: RHEA ID (without 'RHEA:' prefix)

        Returns:
            EQUATION string or None if not found
        """
        return self._equations.get(rhea_id)

    def get_chebi_ids(self, rhea_id: str) -> Tuple[List[str], List[str]]:
        """Get parsed ChEBI IDs for a RHEA ID.

        Args:
            rhea_id: RHEA ID (without 'RHEA:' prefix)

        Returns:
            Tuple of (left_chebi_ids, right_chebi_ids)

        Examples:
            >>> lookup = EquationChEBILookup()
            >>> lookup._equations["10000"] = "CHEBI:16459 + CHEBI:15377 = CHEBI:31011"
            >>> lookup.get_chebi_ids("10000")
            (['CHEBI:16459', 'CHEBI:15377'], ['CHEBI:31011'])
            >>> lookup.get_chebi_ids("99999")
            ([], [])
        """
        # Check cache first
        if rhea_id in self._parsed_cache:
            return self._parsed_cache[rhea_id]

        equation = self._equations.get(rhea_id)
        if not equation:
            return [], []

        result = parse_equation_chebi_ids(equation)
        self._parsed_cache[rhea_id] = result
        return result

    def __len__(self) -> int:
        """Return number of loaded equations."""
        return len(self._equations)
