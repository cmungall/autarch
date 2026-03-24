"""Shared biochemical-context helpers for GO wrapper classifiers."""

from __future__ import annotations

from autarch.datamodel import Participant, PolymerType
from autarch.moiety import Moiety
from autarch.ontology.lipid_utils import carbon_count, nitrogen_count, oxygen_count, phosphate_count

PROTEIN_POLYMERS = {PolymerType.PROTEIN, PolymerType.PEPTIDE, PolymerType.POLYPEPTIDE}


def is_protein_like(participant: Participant) -> bool:
    """Return True for explicit or structure-implied protein participants."""
    if participant.polymer_type in PROTEIN_POLYMERS:
        return True
    if participant.smiles is None or participant.smiles.count("*") < 2:
        return False
    return carbon_count(participant) >= 2 and nitrogen_count(participant) >= 1


def is_steroid_like(participant: Participant) -> bool:
    """Return True for compact polycyclic steroid-like scaffolds."""
    mol = participant.get_mol()
    if mol is None or participant.is_thioester():
        return False
    if phosphate_count(participant) > 0:
        return False
    carbons = carbon_count(participant)
    if carbons < 17 or carbons > 40:
        return False
    return mol.GetRingInfo().NumRings() >= 4


def is_amino_acid_like(participant: Participant) -> bool:
    """Return True for free amino-acid-like small molecules."""
    if participant.is_thioester() or is_protein_like(participant):
        return False
    if not participant.has_moiety(Moiety.AMINO_ACID):
        return False
    return carbon_count(participant) >= 2 and nitrogen_count(participant) >= 1 and oxygen_count(participant) >= 2
