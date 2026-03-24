"""Parser for RHEA reaction labels to extract stoichiometry information.

This module parses human-readable RHEA reaction equations/labels to extract:
- Variable stoichiometry (n, 2n, n+1) for participants
- Polymer indices for polymeric species (n, n+1, n-1)
- Polymer types (RNA, glucan, chitin, etc.)
- Linkage information (1->4, 2->8, etc.)

Examples of RHEA labels with polymer notation:
    - "n ATP + n H2O = n ADP + n phosphate + n H(+)"
    - "[(1->4)-glucosyl](n) + n H2O = [(1->4)-glucosyl](n-1) + glucose"
    - "tRNA(n) + ATP = tRNA(n+1) + diphosphate"
"""

import re
from dataclasses import dataclass
from typing import Optional


@dataclass
class ParsedParticipant:
    """A participant parsed from a reaction label.

    Attributes:
        name: The molecule name (may include polymer notation)
        stoichiometry: Variable stoichiometry like "n", "2n", or None for count=1
        polymer_index: For polymers, the index like "n", "n+1", "n-1"
        base_name: Name without the polymer index suffix
        polymer_type: Inferred polymer type (e.g., "rna", "glucan")
        linkage: Glycosidic/phosphodiester linkage (e.g., "1->4")
        monomer: Description of the repeating unit
        location: Cellular location like "in", "out", "cytoplasm" (from transport reactions)
    """
    name: str
    stoichiometry: Optional[str] = None
    polymer_index: Optional[str] = None
    base_name: Optional[str] = None
    polymer_type: Optional[str] = None
    linkage: Optional[str] = None
    monomer: Optional[str] = None
    location: Optional[str] = None

    def __post_init__(self):
        if self.base_name is None:
            self.base_name = self.name
        # Infer polymer type from name (even without polymer_index for tRNA, proteins, etc.)
        if not self.polymer_type:
            self.polymer_type, self.linkage, self.monomer = infer_polymer_info(self.base_name or self.name)


def infer_polymer_info(name: str) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Infer polymer type, linkage, and monomer from a molecule name.

    Args:
        name: Molecule name, possibly with bracket notation

    Returns:
        Tuple of (polymer_type, linkage, monomer)

    Examples:
        >>> infer_polymer_info("RNA")
        ('rna', None, 'ribonucleotide')
        >>> infer_polymer_info("tRNA")
        ('trna', None, 'ribonucleotide')
        >>> infer_polymer_info("[(1->4)-alpha-D-glucosyl]")
        ('glucan', '1->4', 'alpha-D-glucosyl')
        >>> infer_polymer_info("[(1->4)-N-acetyl-beta-D-glucosaminyl]")
        ('chitin', '1->4', 'N-acetyl-beta-D-glucosaminyl')
        >>> infer_polymer_info("polysulfur")
        ('polysulfur', None, 'sulfur')
        >>> infer_polymer_info("ATP")
        (None, None, None)
    """
    name_lower = name.lower()

    # Nucleic acids
    if name_lower == 'rna' or name_lower.endswith('-rna'):
        return 'rna', None, 'ribonucleotide'
    if name_lower == 'dna' or name_lower.endswith('-dna'):
        return 'dna', None, 'deoxyribonucleotide'
    if name_lower == 'trna' or 'transfer rna' in name_lower or name_lower.startswith('trna('):
        return 'trna', None, 'ribonucleotide'
    # Aminoacyl-tRNA patterns like "L-tyrosyl-tRNA(Tyr)"
    if '-trna(' in name_lower or '-trna' == name_lower[-5:]:
        return 'trna', None, 'ribonucleotide'
    if name_lower == 'mrna' or 'messenger rna' in name_lower:
        return 'mrna', None, 'ribonucleotide'

    # Proteins
    if 'protein' in name_lower or 'peptide' in name_lower:
        return 'protein', None, 'amino acid'

    # Polysulfur
    if 'polysulfur' in name_lower or name_lower == 'sulfur':
        return 'polysulfur', None, 'sulfur'

    # Bracketed polymers - extract linkage and monomer
    bracket_match = re.match(r'\[(.+)\]', name)
    if bracket_match:
        content = bracket_match.group(1)

        # Extract linkage like (1->4), (2->8), (2->6)
        linkage_match = re.search(r'\((\d+)->(\d+)\)', content)
        linkage = f"{linkage_match.group(1)}->{linkage_match.group(2)}" if linkage_match else None

        # Remove linkage from content to get monomer
        monomer = re.sub(r'\(\d+->\d+\)-?', '', content).strip()

        # Infer polymer type from monomer
        monomer_lower = monomer.lower()

        # Glucans (glucose polymers)
        if 'glucosyl' in monomer_lower and 'n-acetyl' not in monomer_lower:
            return 'glucan', linkage, monomer

        # Chitin (N-acetyl-D-glucosamine)
        if 'n-acetyl' in monomer_lower and 'glucosaminyl' in monomer_lower:
            return 'chitin', linkage, monomer

        # Peptidoglycan (GlcNAc-MurNAc)
        if 'glcnac' in monomer_lower and 'mur' in monomer_lower:
            return 'peptidoglycan', linkage, monomer

        # Sialic acid polymers (neuraminosyl)
        if 'neuraminosyl' in monomer_lower or 'sialyl' in monomer_lower:
            return 'sialic_acid_polymer', linkage, monomer

        # Galacturonan (pectin)
        if 'galacturonosyl' in monomer_lower:
            return 'galacturonan', linkage, monomer

        # Fructan
        if 'fructo' in monomer_lower:
            return 'fructan', linkage, monomer

        # Phosphate chains
        if 'phospho' in monomer_lower or 'phosphate' in monomer_lower:
            if 'glyceryl' in monomer_lower or 'ribityl' in monomer_lower:
                return 'teichoic_acid', linkage, monomer
            return 'phosphate_chain', linkage, monomer

        # Generic polysaccharide fallback
        if any(sugar in monomer_lower for sugar in ['osyl', 'ose', 'uronyl']):
            return 'polysaccharide', linkage, monomer

        return 'other', linkage, monomer

    return None, None, None


def parse_stoichiometry_from_term(term: str) -> tuple[Optional[str], str]:
    """Parse stoichiometry prefix from a reaction term.

    Args:
        term: A term like "n ATP", "2n H(+)", "(n-1) sulfur", "2 H2O", or just "ATP"

    Returns:
        Tuple of (stoichiometry, molecule_name)
        stoichiometry is None for count=1, "n" for n, "2n" for 2n, etc.

    Examples:
        >>> parse_stoichiometry_from_term("n ATP")
        ('n', 'ATP')
        >>> parse_stoichiometry_from_term("2n H(+)")
        ('2n', 'H(+)')
        >>> parse_stoichiometry_from_term("(n-1) sulfur")
        ('n-1', 'sulfur')
        >>> parse_stoichiometry_from_term("(n+1) product")
        ('n+1', 'product')
        >>> parse_stoichiometry_from_term("2 H2O")
        (None, '2 H2O')
        >>> parse_stoichiometry_from_term("ATP")
        (None, 'ATP')
    """
    term = term.strip()

    # Match parenthesized stoichiometry: "(n-1)", "(n+1)", "(2n)"
    # Pattern: (digit*n+/-digit?) followed by space and molecule
    match = re.match(r'^\((\d*n(?:[+\-]\d+)?)\)\s+(.+)$', term)
    if match:
        return match.group(1), match.group(2).strip()

    # Match bare stoichiometry: "n", "2n", "n+1", etc.
    # Pattern: optional digit + 'n' + optional '+/-' + optional digit
    match = re.match(r'^(\d*n(?:[+\-]\d+)?)\s+(.+)$', term)
    if match:
        return match.group(1), match.group(2).strip()

    return None, term


def parse_polymer_index(name: str) -> tuple[str, Optional[str]]:
    """Extract polymer index from a molecule name.

    Args:
        name: Molecule name potentially with polymer index like "tRNA(n+1)"

    Returns:
        Tuple of (base_name, polymer_index)

    Examples:
        >>> parse_polymer_index("tRNA(n)")
        ('tRNA', 'n')
        >>> parse_polymer_index("tRNA(n+1)")
        ('tRNA', 'n+1')
        >>> parse_polymer_index("[(1->4)-glucosyl](n)")
        ('[(1->4)-glucosyl]', 'n')
        >>> parse_polymer_index("[(1->4)-glucosyl](n-1)")
        ('[(1->4)-glucosyl]', 'n-1')
        >>> parse_polymer_index("ATP")
        ('ATP', None)
        >>> parse_polymer_index("H(+)")
        ('H(+)', None)
    """
    # Handle bracketed names like [(1->4)-glucosyl](n)
    # The bracket content may contain parentheses, so match the closing ] first
    bracket_match = re.match(r'^(\[[^\]]+\])\((n(?:[+\-]\d+)?)\)$', name)
    if bracket_match:
        return bracket_match.group(1), bracket_match.group(2)

    # Pattern for simple polymer index at end: (n), (n+1), (n-1), (n+2), etc.
    # But NOT things like H(+) or Ca(2+)
    # Polymer indices are specifically: n, n+digit, n-digit
    match = re.match(r'^([A-Za-z][A-Za-z0-9\-]*)\((n(?:[+\-]\d+)?)\)$', name)
    if match:
        return match.group(1), match.group(2)

    return name, None


def split_on_reaction_arrow(label: str) -> Optional[tuple[str, str]]:
    """Split a reaction label on the reaction arrow, ignoring arrows inside brackets.

    Args:
        label: Reaction equation string

    Returns:
        Tuple of (left_side, right_side) or None if no arrow found

    Examples:
        >>> split_on_reaction_arrow("A + B = C + D")
        ('A + B ', ' C + D')
        >>> split_on_reaction_arrow("[(1->4)-glucosyl](n) = [(1->4)-glucosyl](n-1)")
        ('[(1->4)-glucosyl](n) ', ' [(1->4)-glucosyl](n-1)')
        >>> split_on_reaction_arrow("A => B")
        ('A ', ' B')
        >>> split_on_reaction_arrow("A <=> B")
        ('A ', ' B')
    """
    # Find reaction arrows outside of brackets
    # Arrows to look for: <=>, =>, <=, ->, =
    arrows = ['<=>', '=>', '<=', '=']  # Order matters - check longer first

    bracket_depth = 0
    paren_depth = 0
    i = 0

    while i < len(label):
        char = label[i]

        if char == '[':
            bracket_depth += 1
        elif char == ']':
            bracket_depth -= 1
        elif char == '(':
            paren_depth += 1
        elif char == ')':
            paren_depth -= 1
        elif bracket_depth == 0 and paren_depth == 0:
            # Check for arrows at this position
            for arrow in arrows:
                if label[i:i+len(arrow)] == arrow:
                    # Don't match -> as it's used in linkage notation like (1->4)
                    # But do match it if it's surrounded by spaces
                    if arrow == '=':
                        # Simple = is always a valid separator
                        return label[:i], label[i+len(arrow):]
                    elif arrow in ['=>', '<=', '<=>']:
                        return label[:i], label[i+len(arrow):]

        i += 1

    return None


def smart_split_on_plus(text: str) -> list[str]:
    """Split text on ' + ' while preserving content inside brackets.

    Args:
        text: String to split

    Returns:
        List of terms

    Examples:
        >>> smart_split_on_plus("A + B + C")
        ['A', 'B', 'C']
        >>> smart_split_on_plus("[(1->4)-glucosyl](n) + H2O")
        ['[(1->4)-glucosyl](n)', 'H2O']
        >>> smart_split_on_plus("n ATP + n H2O")
        ['n ATP', 'n H2O']
    """
    terms = []
    current_term = []
    bracket_depth = 0
    paren_depth = 0
    i = 0

    while i < len(text):
        char = text[i]

        if char == '[':
            bracket_depth += 1
            current_term.append(char)
        elif char == ']':
            bracket_depth -= 1
            current_term.append(char)
        elif char == '(':
            paren_depth += 1
            current_term.append(char)
        elif char == ')':
            paren_depth -= 1
            current_term.append(char)
        elif char == '+' and bracket_depth == 0 and paren_depth == 0:
            # Check if this is " + " (with spaces)
            if i > 0 and i < len(text) - 1:
                before = text[i-1] if i > 0 else ''
                after = text[i+1] if i < len(text) - 1 else ''
                if before == ' ' and after == ' ':
                    # This is a separator
                    term = ''.join(current_term).strip()
                    if term:
                        terms.append(term)
                    current_term = []
                    i += 1  # Skip the space after +
                else:
                    current_term.append(char)
            else:
                current_term.append(char)
        else:
            current_term.append(char)

        i += 1

    # Add the last term
    term = ''.join(current_term).strip()
    if term:
        terms.append(term)

    return terms


def parse_location_annotation(name: str) -> tuple[str, Optional[str]]:
    """Extract location annotation from a molecule name.

    Location annotations like (in), (out), (cytoplasm) appear at the
    end of molecule names in transport reactions.

    Args:
        name: Molecule name that may have location annotation

    Returns:
        Tuple of (clean_name, location) where location is None if not present

    Examples:
        >>> parse_location_annotation("sulfate(out)")
        ('sulfate', 'out')
        >>> parse_location_annotation("Mg(2+)(in)")
        ('Mg(2+)', 'in')
        >>> parse_location_annotation("ATP")
        ('ATP', None)
        >>> parse_location_annotation("calcium(2+)(cytoplasm)")
        ('calcium(2+)', 'cytoplasm')
    """
    # Common location annotations in RHEA
    location_patterns = [
        r'\(in\)$',
        r'\(out\)$',
        r'\(cytoplasm\)$',
        r'\(periplasm\)$',
        r'\(extracellular\)$',
        r'\(matrix\)$',
        r'\(intermembrane space\)$',
    ]

    for pattern in location_patterns:
        match = re.search(pattern, name, re.IGNORECASE)
        if match:
            location = match.group(0)[1:-1]  # Remove parentheses
            clean_name = name[:match.start()]
            return clean_name, location.lower()

    return name, None


def parse_reaction_label(label: str) -> tuple[list[ParsedParticipant], list[ParsedParticipant]]:
    """Parse a RHEA reaction label into left and right participants.

    Args:
        label: A RHEA reaction equation like:
            "n ATP + n H2O = n ADP + n phosphate"
            "tRNA(n) + ATP = tRNA(n+1) + diphosphate"

    Returns:
        Tuple of (left_participants, right_participants)

    Examples:
        >>> left, right = parse_reaction_label("n ATP + n H2O = n ADP + n Pi")
        >>> len(left)
        2
        >>> left[0].name
        'ATP'
        >>> left[0].stoichiometry
        'n'
        >>> right[0].name
        'ADP'
        >>> right[0].stoichiometry
        'n'

        >>> left, right = parse_reaction_label("tRNA(n) + ATP = tRNA(n+1) + PPi")
        >>> left[0].name
        'tRNA'
        >>> left[0].polymer_index
        'n'
        >>> right[0].name
        'tRNA'
        >>> right[0].polymer_index
        'n+1'
    """
    # Split on reaction arrow - must find the arrow outside of brackets
    parts = split_on_reaction_arrow(label)
    if parts is None:
        return [], []

    left_str, right_str = parts

    def parse_side(side_str: str) -> list[ParsedParticipant]:
        """Parse one side of the reaction equation."""
        participants = []

        # Split by '+' but preserve content inside brackets
        # We need to be careful not to split on + inside brackets like [(1->4)]
        terms = smart_split_on_plus(side_str.strip())

        for term in terms:
            term = term.strip()
            if not term:
                continue

            # Extract stoichiometry
            stoich, name = parse_stoichiometry_from_term(term)

            # Extract location annotation (e.g., sulfate(out) -> sulfate, "out")
            name_without_loc, location = parse_location_annotation(name)

            # Extract polymer index
            base_name, polymer_idx = parse_polymer_index(name_without_loc)

            participants.append(ParsedParticipant(
                name=base_name,
                stoichiometry=stoich,
                polymer_index=polymer_idx,
                base_name=base_name,
                location=location
            ))

        return participants

    return parse_side(left_str), parse_side(right_str)


def has_polymer_notation(label: str) -> bool:
    """Check if a label contains polymer (n) notation.

    Args:
        label: RHEA reaction label

    Returns:
        True if label contains (n), (n+1), (n-1), or "n " stoichiometry

    Examples:
        >>> has_polymer_notation("n ATP + n H2O = n ADP + n Pi")
        True
        >>> has_polymer_notation("tRNA(n) + ATP = tRNA(n+1) + PPi")
        True
        >>> has_polymer_notation("(n-1) sulfur + H2S = n sulfur + H2")
        True
        >>> has_polymer_notation("ATP + H2O = ADP + Pi")
        False
    """
    # Check for polymer indices like (n), (n+1), (n-1) at end of molecule names
    # e.g., tRNA(n), polymer(n+1)
    if re.search(r'\w\(n(?:[+\-]\d+)?\)', label):
        return True

    # Check for stoichiometry coefficients like "(n-1) sulfur", "(n+1) product"
    if re.search(r'\(\d*n(?:[+\-]\d+)?\)\s+\w', label):
        return True

    # Check for variable stoichiometry like "n ATP", "2n H(+)"
    if re.search(r'\b\d*n\s+[A-Za-z\[\(]', label):
        return True

    return False


def extract_polymer_info(label: str) -> dict[str, object]:
    """Extract detailed polymer information from a label.

    Args:
        label: RHEA reaction label

    Returns:
        Dictionary with:
            - has_polymer: bool
            - polymer_participants: list of names with polymer indices
            - variable_stoich_participants: list of names with n-stoichiometry

    Examples:
        >>> info = extract_polymer_info("tRNA(n) + ATP = tRNA(n+1) + PPi")
        >>> info['has_polymer']
        True
        >>> info['polymer_participants']
        [('tRNA', 'n'), ('tRNA', 'n+1')]
    """
    polymer_participants: list[tuple[str, str]] = []
    variable_stoich_participants: list[tuple[str, str]] = []

    if not has_polymer_notation(label):
        return {
            "has_polymer": False,
            "polymer_participants": polymer_participants,
            "variable_stoich_participants": variable_stoich_participants,
        }

    left, right = parse_reaction_label(label)

    for p in left + right:
        if p.polymer_index:
            polymer_participants.append((p.base_name or p.name, p.polymer_index))
        if p.stoichiometry:
            variable_stoich_participants.append((p.name, p.stoichiometry))

    return {
        "has_polymer": True,
        "polymer_participants": polymer_participants,
        "variable_stoich_participants": variable_stoich_participants,
    }
