"""Utilities for transport-classifier substrate family detection.

These helpers keep transport subclasses focused on the transported substrate
family rather than repeating transport-pair bookkeeping in each classifier.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.moiety import Moiety
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NH4,
    CHEBI_PHOSPHATE,
)

CHEBI_POLYPHOSPHATE = "CHEBI:16838"
CHEBI_SODIUM = "CHEBI:29101"
CHEBI_POTASSIUM = "CHEBI:29103"
KNOWN_INORGANIC_CATION_IDS = {
    CHEBI_SODIUM,
    CHEBI_POTASSIUM,
    "CHEBI:29105",  # lithium(1+)
    "CHEBI:29106",  # cesium(1+)
    "CHEBI:29107",  # rubidium(1+)
    "CHEBI:29108",  # calcium(2+)
    "CHEBI:29033",  # Fe(2+)
    "CHEBI:29034",  # Fe(3+)
    "CHEBI:29104",  # magnesium(2+)
}
KNOWN_TRANSITION_METAL_IDS = {
    "CHEBI:29033",  # Fe(2+)
    "CHEBI:29034",  # Fe(3+)
    "CHEBI:29108",  # calcium(2+)
    "CHEBI:29104",  # magnesium(2+)
}
CHEBI_GENERIC_AMINO_ACID_PROTEIN = "CHEBI:83228"
CHEBI_GENERIC_DIPEPTIDE = "CHEBI:90799"
CHEBI_GENERIC_GLUTATHIONE_CHAIN = "CHEBI:131728"
CHEBI_GENERIC_QUATERNARY_AMMONIUM_CARRIER = "CHEBI:35267"
CHEBI_GENERIC_ACYL_COA_CARRIER = "CHEBI:58342"
CHEBI_GENERIC_COBALAMIN_CARRIER = "CHEBI:140785"
CHEBI_GENERIC_GLUTATHIONE_CARRIER = "CHEBI:90779"

TRANSPORT_SPECTATOR_IDS = {
    CHEBI_ATP,
    CHEBI_ADP,
    CHEBI_GTP,
    CHEBI_GDP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    CHEBI_POLYPHOSPHATE,
}

GENERIC_AMINO_PEPTIDE_IDS = {
    CHEBI_GENERIC_AMINO_ACID_PROTEIN,
    CHEBI_GENERIC_DIPEPTIDE,
    CHEBI_GENERIC_GLUTATHIONE_CHAIN,
}

SUPPORTED_POLAR_AMINO_ACID_IDS = {
    "CHEBI:62031",
    "CHEBI:32682",
    "CHEBI:46911",
    "CHEBI:32551",
}

GENERIC_TRANSPORT_CARRIER_PLACEHOLDER_IDS = {
    CHEBI_GENERIC_QUATERNARY_AMMONIUM_CARRIER,
    CHEBI_GENERIC_ACYL_COA_CARRIER,
    CHEBI_GENERIC_COBALAMIN_CARRIER,
    CHEBI_GENERIC_GLUTATHIONE_CARRIER,
}

CARBOHYDRATE_POLYMER_TYPES = {
    PolymerType.GLUCAN,
    PolymerType.STARCH,
    PolymerType.GLYCOGEN,
    PolymerType.CELLULOSE,
    PolymerType.CHITIN,
    PolymerType.CHITOSAN,
    PolymerType.PEPTIDOGLYCAN,
    PolymerType.SIALIC_ACID_POLYMER,
    PolymerType.GALACTURONAN,
    PolymerType.FRUCTAN,
    PolymerType.POLYSACCHARIDE,
    PolymerType.TEICHOIC_ACID,
}

AMINO_PEPTIDE_POLYMER_TYPES = {
    PolymerType.PROTEIN,
    PolymerType.PEPTIDE,
    PolymerType.POLYPEPTIDE,
}

AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
RING_OXYGEN_PATTERN = Chem.MolFromSmarts("[O;R]")
ALPHA_AMINO_ACID_PATTERN = Chem.MolFromSmarts(
    "[NX3,NX4+;!$([N+](C)(C)C)][CX4][CX3](=[OX1])[O-,$([OX2H1])]"
)
CARBOXYL_PATTERN = Chem.MolFromSmarts("[C;X3](=[O;X1])[O;H1,X1-]")


def transported_pairs(reaction: Reaction) -> list[tuple[Participant, Participant]]:
    """Return transported participant pairs across the membrane."""
    pairs: list[tuple[Participant, Participant]] = []
    for left_participant in reaction.left_participants:
        for right_participant in reaction.right_participants:
            if (
                left_participant.is_same_molecule(right_participant)
                and left_participant.location != right_participant.location
            ):
                pairs.append((left_participant, right_participant))
    return pairs


def reactive_transported_pairs(reaction: Reaction) -> list[tuple[Participant, Participant]]:
    """Return transported pairs excluding ATPase-coupling spectators."""
    return [
        (left_participant, right_participant)
        for left_participant, right_participant in transported_pairs(reaction)
        if left_participant.chebi_id not in TRANSPORT_SPECTATOR_IDS
        and left_participant.chebi_id not in GENERIC_TRANSPORT_CARRIER_PLACEHOLDER_IDS
    ]


def coupled_transported_pairs(reaction: Reaction) -> list[tuple[Participant, Participant]]:
    """Return transported pairs, retaining H+ for coupled transport logic."""
    retained_spectators = TRANSPORT_SPECTATOR_IDS - {CHEBI_H_PLUS}
    return [
        (left_participant, right_participant)
        for left_participant, right_participant in transported_pairs(reaction)
        if left_participant.chebi_id not in retained_spectators
        and left_participant.chebi_id not in GENERIC_TRANSPORT_CARRIER_PLACEHOLDER_IDS
    ]


def atom_count(participant: Participant, atomic_number: int) -> int:
    """Count atoms of one element in a participant."""
    mol = participant.get_mol()
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_number)


def formal_charge(participant: Participant) -> int:
    """Return the total formal charge of a participant."""
    mol = participant.get_mol()
    if mol is None:
        return 0
    return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())


def is_inorganic_anion(participant: Participant) -> bool:
    """Return True for directly represented inorganic anions."""
    if participant.chebi_id == CHEBI_POLYPHOSPHATE:
        return True
    mol = participant.get_mol()
    if mol is None:
        return False
    return atom_count(participant, 6) == 0 and formal_charge(participant) < 0


def is_inorganic_cation(participant: Participant) -> bool:
    """Return True for directly represented inorganic cations."""
    if participant.chebi_id in KNOWN_INORGANIC_CATION_IDS:
        return True

    mol = participant.get_mol()
    if mol is None:
        return False

    carbon_count = atom_count(participant, 6)
    charge = formal_charge(participant)
    return carbon_count == 0 and charge > 0


def is_metal_ion(participant: Participant) -> bool:
    """Return True for transported metal ions, excluding hydron and ammonium."""
    return is_inorganic_cation(participant) and participant.chebi_id not in {CHEBI_H_PLUS, CHEBI_NH4}


def is_sodium_ion(participant: Participant) -> bool:
    """Return True for sodium ion."""
    return participant.chebi_id == CHEBI_SODIUM


def is_transition_metal_ion(participant: Participant) -> bool:
    """Return True for the supported transition-metal-like transported ions."""
    return participant.chebi_id in KNOWN_TRANSITION_METAL_IDS


def is_neutral_l_amino_acid(participant: Participant) -> bool:
    """Return True for neutral alpha-amino-acid transporter substrates."""
    mol = participant.get_mol()
    if mol is None or ALPHA_AMINO_ACID_PATTERN is None:
        return False
    if not mol.HasSubstructMatch(ALPHA_AMINO_ACID_PATTERN):
        return False
    if any(atom.GetAtomicNum() not in {1, 6, 7, 8, 16} for atom in mol.GetAtoms()):
        return False
    if _count_pattern(participant, CARBOXYL_PATTERN) != 1:
        return False
    positively_charged_nitrogens = sum(
        1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() > 0
    )
    return positively_charged_nitrogens <= 1


def is_amino_acid_or_peptide(participant: Participant) -> bool:
    """Return True for amino-acid and peptide-like transported substrates."""
    if participant.polymer_type in AMINO_PEPTIDE_POLYMER_TYPES:
        return True
    if participant.chebi_id in GENERIC_AMINO_PEPTIDE_IDS:
        return True

    mol = participant.get_mol()
    if mol is None:
        return False

    if any(atom.GetAtomicNum() not in {1, 6, 7, 8, 16} for atom in mol.GetAtoms()):
        return False

    carbon_count = atom_count(participant, 6)
    if carbon_count < 2:
        return False

    if ALPHA_AMINO_ACID_PATTERN is not None and mol.HasSubstructMatch(ALPHA_AMINO_ACID_PATTERN):
        return True

    if AMIDE_PATTERN is None:
        return False

    amide_count = len(mol.GetSubstructMatches(AMIDE_PATTERN))
    return amide_count > 0 and atom_count(participant, 7) >= 1


def is_supported_polar_amino_acid(participant: Participant) -> bool:
    """Return True for the supported polar amino-acid transporter substrates."""
    return participant.chebi_id in SUPPORTED_POLAR_AMINO_ACID_IDS


def is_carbohydrate_or_derivative(participant: Participant) -> bool:
    """Return True for carbohydrate-like transported substrates."""
    if participant.polymer_type in CARBOHYDRATE_POLYMER_TYPES:
        return True

    mol = participant.get_mol()
    if mol is None:
        return False

    carbon_count = atom_count(participant, 6)
    if carbon_count < 3:
        return False

    if participant.is_glycoside():
        return True

    if RING_OXYGEN_PATTERN is None:
        return False

    ring_oxygen_count = len(mol.GetSubstructMatches(RING_OXYGEN_PATTERN))
    hydroxyl_like = participant.has_moiety(Moiety.HYDROXYL)
    return ring_oxygen_count > 0 and hydroxyl_like and atom_count(participant, 7) == 0


def is_carboxylic_acid_or_derivative(participant: Participant) -> bool:
    """Return True for transported carboxylic acids or carboxylates."""
    mol = participant.get_mol()
    if mol is None:
        return False
    return participant.has_moiety(Moiety.CARBOXYL)


def is_dicarboxylic_acid(participant: Participant) -> bool:
    """Return True for carboxylic acids with at least two carboxyl groups."""
    return is_carboxylic_acid_or_derivative(participant) and _count_pattern(participant, CARBOXYL_PATTERN) >= 2


def is_sulfur_compound(participant: Participant) -> bool:
    """Return True for directly represented sulfur-containing compounds."""
    if participant.chebi_id == CHEBI_H_PLUS:
        return False
    return atom_count(participant, 16) > 0


def _count_pattern(participant: Participant, pattern: Chem.Mol | None) -> int:
    """Count SMARTS matches in one participant."""
    mol = participant.get_mol()
    if mol is None or pattern is None:
        return 0
    return len(mol.GetSubstructMatches(pattern))
