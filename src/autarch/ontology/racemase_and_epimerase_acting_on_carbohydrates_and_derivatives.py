"""racemase and epimerase activity, acting on carbohydrates and derivatives.

Catalysis of a reaction that alters the configuration of one or more chiral centers in a carbohydrate or carbohydrate derivative.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.racemase_and_epimerase import RacemaseAndEpimerase
from autarch.stereochemistry import are_stereoisomers


class RacemaseAndEpimeraseActingOnCarbohydratesAndDerivatives(RacemaseAndEpimerase):
    """racemase and epimerase activity, acting on carbohydrates and derivatives.

    Catalysis of a reaction that alters the configuration of one or more chiral
    centers in a carbohydrate or carbohydrate derivative.
    """

    GO_ID = "GO:0016857"
    EC_NUMBER_PREFIX = "5.1.3.-"
    EXCLUDED_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_ATP,
        CHEBI_ADP,
        CHEBI_CO2,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for carbohydrate stereochemical inversion without cofactors."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not carbohydrate racemase chemistry",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External cofactors indicate a different reaction class",
            )

        left_candidates = [
            participant
            for participant in reaction.left_participants
            if self._is_carbohydrate_derivative(participant)
        ]
        right_candidates = [
            participant
            for participant in reaction.right_participants
            if self._is_carbohydrate_derivative(participant)
        ]
        if len(left_candidates) != 1 or len(right_candidates) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Requires one structured carbohydrate-like participant on each side",
            )

        substrate = left_candidates[0]
        product = right_candidates[0]
        if not self._same_connectivity_ignore_stereo(substrate, product):
            return ClassificationResult(
                is_member=False,
                explanation="Carbohydrate epimerization conserves molecular connectivity",
            )

        if not self._has_stereochemical_inversion(substrate, product):
            return ClassificationResult(
                is_member=False,
                explanation="No stereochemical inversion detected for the carbohydrate-like scaffold",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Carbohydrate racemase/epimerase: stereochemical inversion of a carbohydrate-like scaffold",
        )

    @staticmethod
    def _is_carbohydrate_derivative(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
        sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
        if carbon_count < 4 or oxygen_count < 2 or sulfur_count > 0:
            return False
        if not Chem.FindMolChiralCenters(mol, includeUnassigned=True):
            return False
        has_ring_oxygen = any(atom.GetAtomicNum() == 8 and atom.IsInRing() for atom in mol.GetAtoms())
        oxygen_rich = oxygen_count * 2 >= carbon_count
        return has_ring_oxygen or phosphorus_count > 0 or oxygen_rich

    @staticmethod
    def _same_connectivity_ignore_stereo(left: Participant, right: Participant) -> bool:
        if left.inchi and right.inchi:
            return (
                left.inchi.split("/t")[0].split("/m")[0].split("/s")[0]
                == right.inchi.split("/t")[0].split("/m")[0].split("/s")[0]
            )
        left_mol = left.get_mol()
        right_mol = right.get_mol()
        if left_mol is None or right_mol is None:
            return False
        left_copy = Chem.Mol(left_mol)
        right_copy = Chem.Mol(right_mol)
        Chem.RemoveStereochemistry(left_copy)
        Chem.RemoveStereochemistry(right_copy)
        return Chem.MolToSmiles(left_copy, isomericSmiles=False) == Chem.MolToSmiles(
            right_copy, isomericSmiles=False
        )

    @staticmethod
    def _has_stereochemical_inversion(left: Participant, right: Participant) -> bool:
        if left.inchi and right.inchi:
            is_stereo, _ = are_stereoisomers(left.inchi, right.inchi)
            if is_stereo:
                return True
        left_mol = left.get_mol()
        right_mol = right.get_mol()
        if left_mol is None or right_mol is None:
            return False
        left_iso = Chem.MolToSmiles(left_mol, isomericSmiles=True)
        right_iso = Chem.MolToSmiles(right_mol, isomericSmiles=True)
        left_copy = Chem.Mol(left_mol)
        right_copy = Chem.Mol(right_mol)
        Chem.RemoveStereochemistry(left_copy)
        Chem.RemoveStereochemistry(right_copy)
        left_flat = Chem.MolToSmiles(left_copy, isomericSmiles=False)
        right_flat = Chem.MolToSmiles(right_copy, isomericSmiles=False)
        return left_flat == right_flat and left_iso != right_iso
