"""Catalysis of the hydrolysis of any acid halide bond in substances containing halogen atoms in organic linkage.

This EC 3.8.1 branch covers hydrolytic dehalogenation of carbon-halogen bonds in
organic substrates. The defining signature is water-dependent loss of a C-halide
bond with release of a halide ion and formation of a more oxygenated organic
product.
"""

from __future__ import annotations

from abc import abstractmethod

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.hydrolase import Hydrolase

HALIDE_CHEBI_IDS = {
    "CHEBI:17996",  # chloride
    "CHEBI:24061",  # fluoride
    "CHEBI:29227",  # bromide
    "CHEBI:16788",  # iodide
}


class HydrolaseActingOnAcidHalideBondsInCHalideCompounds(Hydrolase):
    """Legacy implementation for hydrolase acting on halide bonds in C-halide compounds.

    Catalysis of the hydrolysis of any acid halide bond in substances
    containing halogen atoms in organic linkage.

    The class is defined by water-dependent cleavage of an organic carbon-halogen
    bond, release of a halide ion, and loss of carbon-halogen connectivity in the
    substantive substrate scaffold.
    """

    EC_NUMBER_PREFIX = "3.8.1.-"

    @abstractmethod
    def _implementation_only(self) -> None:
        """Mark this legacy implementation class as abstract."""

    EXCLUDED_CHEBIS = {
        CHEBI_ATP,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = HALIDE_CHEBI_IDS | {CHEBI_H_PLUS}
    CARBON_HALIDE_PATTERN = Chem.MolFromSmarts("[c,C][F,Cl,Br,I]")
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1][#6]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolytic dehalogenation of carbon-halogen bonds."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not carbon-halide hydrolysis",
            )

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External redox or ATP coupling indicates another reaction class",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Water is required for hydrolytic dehalogenation",
            )

        if not any(participant.chebi_id in HALIDE_CHEBI_IDS for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No halide ion product detected",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No chemically resolved organic scaffold remains after removing solvent and halide spectators",
            )

        left_c_hal = sum(self._substructure_count(participant, self.CARBON_HALIDE_PATTERN) for participant in left_reactive)
        right_c_hal = sum(self._substructure_count(participant, self.CARBON_HALIDE_PATTERN) for participant in right_reactive)
        if left_c_hal <= right_c_hal or left_c_hal == 0:
            return ClassificationResult(
                is_member=False,
                explanation="No carbon-halogen bond is hydrolytically cleaved",
            )

        left_hydroxyl = sum(
            self._substructure_count(participant, self.HYDROXYL_PATTERN)
            for participant in left_reactive
        )
        right_hydroxyl = sum(
            self._substructure_count(participant, self.HYDROXYL_PATTERN)
            for participant in right_reactive
        )
        if right_hydroxyl <= left_hydroxyl:
            return ClassificationResult(
                is_member=False,
                explanation="The organic products do not show hydroxyl substitution after halide release",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Hydrolase acting on acid halide bonds in C-halide compounds: water replaces an organic halide substituent and releases a halide ion",
        )

    @staticmethod
    def _substructure_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
