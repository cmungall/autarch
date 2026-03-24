"""Shared helpers for lipid and polyphosphate GO classifiers."""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import Participant
from autarch.moiety import Moiety

PHOSPHATE_PATTERN = Chem.MolFromSmarts("[P;X4](=[O;X1])([O;X2,OX1-])([O;X2,OX1-])[O;X2,OX1-]")
DOUBLE_BOND_PATTERN = Chem.MolFromSmarts("[C;X3]=[C;X3]")
RING_OXYGEN_PATTERN = Chem.MolFromSmarts("[O;R]")

GLUCOSE_CHEBIS = {
    "CHEBI:17234",
    "CHEBI:4167",
    "CHEBI:17634",
    "CHEBI:15903",
    "CHEBI:61548",
    "CHEBI:17665",
}


def count_pattern(participant: Participant, pattern: Chem.Mol | None) -> int:
    """Count SMARTS matches in one participant."""
    mol = participant.get_mol()
    if mol is None or pattern is None:
        return 0
    return len(mol.GetSubstructMatches(pattern))


def atom_count(participant: Participant, atomic_number: int) -> int:
    """Count atoms of one element in a participant."""
    mol = participant.get_mol()
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_number)


def carbon_count(participant: Participant) -> int:
    """Count carbons in a participant."""
    return atom_count(participant, 6)


def oxygen_count(participant: Participant) -> int:
    """Count oxygens in a participant."""
    return atom_count(participant, 8)


def nitrogen_count(participant: Participant) -> int:
    """Count nitrogens in a participant."""
    return atom_count(participant, 7)


def sulfur_count(participant: Participant) -> int:
    """Count sulfurs in a participant."""
    return atom_count(participant, 16)


def phosphate_count(participant: Participant) -> int:
    """Count phosphate groups in a participant."""
    return count_pattern(participant, PHOSPHATE_PATTERN)


def double_bond_count(participant: Participant) -> int:
    """Count C=C double bonds in a participant."""
    return count_pattern(participant, DOUBLE_BOND_PATTERN)


def is_glucose_like(participant: Participant) -> bool:
    """Return True for glucose or glucose-phosphate products."""
    if participant.chebi_id in GLUCOSE_CHEBIS:
        return True
    mol = participant.get_mol()
    if mol is None:
        return False
    if carbon_count(participant) != 6 or oxygen_count(participant) < 5:
        return False
    return count_pattern(participant, RING_OXYGEN_PATTERN) >= 1


def is_inositol_phosphate_like(participant: Participant) -> bool:
    """Return True for soluble inositol phosphates."""
    mol = participant.get_mol()
    if mol is None or participant.is_thioester():
        return False
    if carbon_count(participant) != 6 or oxygen_count(participant) < 6:
        return False
    if phosphate_count(participant) < 1:
        return False
    if nitrogen_count(participant) > 0 or sulfur_count(participant) > 0:
        return False
    return mol.GetRingInfo().NumRings() >= 1


def is_glycerophospholipid_like(participant: Participant) -> bool:
    """Return True for glycerophospholipid-like participants."""
    mol = participant.get_mol()
    if mol is None or participant.is_thioester():
        return False
    if carbon_count(participant) < 8 or oxygen_count(participant) < 5:
        return False
    return phosphate_count(participant) >= 1


def is_phosphatidylinositol_like(participant: Participant) -> bool:
    """Return True for phosphatidylinositol and phosphatidylinositol phosphates."""
    mol = participant.get_mol()
    if mol is None or participant.is_thioester():
        return False
    if carbon_count(participant) < 10 or oxygen_count(participant) < 8:
        return False
    return phosphate_count(participant) >= 1


def is_lipid_like(participant: Participant) -> bool:
    """Return True for lipid-like organic participants."""
    if participant.is_thioester():
        return False
    mol = participant.get_mol()
    if mol is None:
        return False
    return carbon_count(participant) >= 8 and (
        participant.is_phosphorylated()
        or participant.has_moiety(Moiety.CARBOXYL)
        or oxygen_count(participant) >= 2
    )


def is_fatty_acyl_coa_like(participant: Participant) -> bool:
    """Return True for fatty acyl-CoA-like thioesters."""
    return participant.is_thioester() and carbon_count(participant) >= 8


def is_prenyl_diphosphate_like(participant: Participant) -> bool:
    """Return True for prenyl diphosphate substrates and products."""
    mol = participant.get_mol()
    if mol is None:
        return False
    if phosphate_count(participant) < 2:
        return False
    if nitrogen_count(participant) > 0 or sulfur_count(participant) > 0:
        return False
    return carbon_count(participant) in {5, 10, 15, 20, 25}
