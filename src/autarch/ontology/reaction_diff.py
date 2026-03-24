"""Ultra-simple reaction analysis based on molecular-level changes.

Key insight: Don't track individual bonds. Instead:
- Hydrolysis = water + fragmentation
- Oxidoreduction = redox cofactors
- Transfer = group moved between molecules
- Etc.
"""

from typing import Dict, List, Set
from collections import defaultdict
from rdkit import Chem

from autarch.datamodel import Reaction, Participant

from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)


class ReactionDiff:
    """
    Ultra-simple reaction analyzer based on molecular-level changes.

    Philosophy: Don't overthink it. Track molecular changes, not bond details.
    """

    def __init__(self, reaction: Reaction):
        self.reaction = reaction

        # Count molecules (not atoms)
        self.n_reactant_molecules = self._count_molecules(reaction.left_participants)
        self.n_product_molecules = self._count_molecules(reaction.right_participants)

        # Fragmentation/fusion
        self.is_fragmentation = self.n_product_molecules > self.n_reactant_molecules
        self.is_fusion = self.n_product_molecules < self.n_reactant_molecules

        # Water detection
        self.has_water_reactant = self._has_water(reaction.left_participants)
        self.has_water_product = self._has_water(reaction.right_participants)

        # Atom types present
        self.reactant_elements = self._get_elements(reaction.left_participants)
        self.product_elements = self._get_elements(reaction.right_participants)

        # Atom counts for balance
        self.reactant_atoms = self._count_atoms(reaction.left_participants)
        self.product_atoms = self._count_atoms(reaction.right_participants)

        # Check for specific cofactors
        self.has_redox_cofactor = self._has_redox_cofactor()
        self.has_atp = self._has_atp()
        self.has_phosphate = "P" in self.reactant_elements

    def _count_molecules(self, participants: List[Participant]) -> int:
        """Count number of molecules, excluding water and H+."""
        skip_chebi = {
            CHEBI_H2O,  # water
            CHEBI_H_PLUS,  # H+
        }

        count = 0
        for stoi in participants:
            if stoi.chebi_id not in skip_chebi:
                count += stoi.count
        return count

    def _has_water(self, participants: List[Participant]) -> bool:
        """Check if water is present."""
        water_chebi = CHEBI_H2O
        for stoi in participants:
            if stoi.chebi_id == water_chebi:
                return True
            mol = stoi.get_mol()
            if mol:
                # Add implicit hydrogens for accurate atom count
                mol_with_h = Chem.AddHs(mol)
                if mol_with_h.GetNumAtoms() == 3:
                    atoms = sorted([a.GetSymbol() for a in mol_with_h.GetAtoms()])
                    if atoms == ["H", "H", "O"]:
                        return True
        return False

    def _get_elements(self, participants: List[Participant]) -> Set[str]:
        """Get all element types present."""
        elements = set()
        for stoi in participants:
            mol = stoi.get_mol()
            if mol:
                # Add implicit hydrogens for accurate element set
                mol_with_h = Chem.AddHs(mol)
                for atom in mol_with_h.GetAtoms():
                    elements.add(atom.GetSymbol())
        return elements

    def _count_atoms(self, participants: List[Participant]) -> Dict[str, int]:
        """Count atoms."""
        counts: Dict[str, int] = defaultdict(int)
        for stoi in participants:
            mol = stoi.get_mol()
            if mol:
                # Add implicit hydrogens for accurate atom count
                mol_with_h = Chem.AddHs(mol)
                for atom in mol_with_h.GetAtoms():
                    counts[atom.GetSymbol()] += stoi.count
        return dict(counts)

    def _has_redox_cofactor(self) -> bool:
        """Check for redox cofactors."""
        redox_cofactors = {
            CHEBI_NAD_PLUS,  # NAD+
            CHEBI_NADH,  # NADH
            CHEBI_NADP_PLUS,  # NADP+
            CHEBI_NADPH,  # NADPH
            CHEBI_FAD,  # FAD
            CHEBI_FADH2,  # FADH2
            CHEBI_O2,  # O2
        }

        for stoi in self.reaction.left_participants + self.reaction.right_participants:
            if stoi.chebi_id in redox_cofactors:
                return True
        return False

    def _has_atp(self) -> bool:
        """Check for ATP/ADP."""
        atp_adp = {
            CHEBI_ATP,  # ATP
            CHEBI_ADP,  # ADP
        }
        for stoi in self.reaction.left_participants:
            if stoi.chebi_id in atp_adp:
                return True
        return False

    def involves_bond_breaking(self, bond_type: str) -> bool:
        """
        Ultra-simple bond breaking detection.
        If water is consumed AND fragmentation occurs AND the atoms are present,
        then this bond type was likely broken.
        """
        if not (self.has_water_reactant and self.is_fragmentation):
            return False

        # Parse bond type to get atoms
        for separator in ["≡", "=", "-", "~"]:
            if separator in bond_type:
                atoms = bond_type.split(separator)
                if len(atoms) == 2:
                    atom1, atom2 = atoms
                    return (
                        atom1 in self.reactant_elements
                        and atom2 in self.reactant_elements
                    )

        return False

    def involves_water_addition(self) -> bool:
        """Water consumed in reaction."""
        return self.has_water_reactant

    def involves_water_removal(self) -> bool:
        """Water produced in reaction."""
        return self.has_water_product

    def get_atom_balance(self) -> Dict[str, int]:
        """Get atom balance."""
        balance = {}
        all_atoms = set(self.reactant_atoms.keys()) | set(self.product_atoms.keys())
        for atom in all_atoms:
            balance[atom] = self.product_atoms.get(atom, 0) - self.reactant_atoms.get(
                atom, 0
            )
        return balance

    def is_balanced(self) -> bool:
        """Check if balanced (allowing for H+ production)."""
        balance = self.get_atom_balance()

        # Allow H imbalance of ±1 (proton production/consumption)
        for atom, change in balance.items():
            if atom == "H" and abs(change) <= 1:
                continue
            if change != 0:
                return False
        return True

    def involves_oxidation_state_change(self) -> bool:
        """
        Simple: if redox cofactors present, it's oxidation.
        """
        return self.has_redox_cofactor

    def is_likely_hydrolase(self) -> bool:
        """
        Ultra-simple hydrolase detection:
        Water consumed + fragmentation + no redox cofactors
        """
        return (
            self.has_water_reactant
            and self.is_fragmentation
            and not self.has_redox_cofactor
        )

    def is_likely_oxidoreductase(self) -> bool:
        """
        Ultra-simple oxidoreductase detection:
        Has redox cofactors
        """
        return self.has_redox_cofactor

    def is_likely_transferase(self) -> bool:
        """
        Ultra-simple transferase detection:
        ATP/ADP present (phosphoryl transfer) or
        No fragmentation but atoms rearranged
        """
        return self.has_atp and not self.is_fragmentation

    def involves_bond_formation(self, bond_type: str) -> bool:
        """
        Simple bond formation detection.
        For compatibility - just returns False since we don't track this anymore.
        """
        return False

    def get_summary(self) -> str:
        """Simple summary."""
        parts = []

        if self.is_fragmentation:
            parts.append(
                f"Fragmentation ({self.n_reactant_molecules}→{self.n_product_molecules})"
            )
        elif self.is_fusion:
            parts.append(
                f"Fusion ({self.n_reactant_molecules}→{self.n_product_molecules})"
            )

        if self.has_water_reactant:
            parts.append("Water consumed")
        if self.has_water_product:
            parts.append("Water produced")

        if self.has_redox_cofactor:
            parts.append("Has redox cofactor")

        if self.has_atp:
            parts.append("Has ATP/ADP")

        if self.has_phosphate:
            parts.append("Contains phosphate")

        return "; ".join(parts) if parts else "No clear pattern"
