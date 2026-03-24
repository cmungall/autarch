"""hydrolase acting on carbon-nitrogen but not peptide bonds in cyclic amides.

Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
cyclic amide.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.hydrolase import Hydrolase

CHEBI_POLYPHOSPHATE = "CHEBI:16838"


class HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmides(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds in cyclic amides.

    Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
    cyclic amide.
    """

    GO_ID = "GO:0016812"
    EC_NUMBER_PREFIX = "3.5.2.-"

    LEFT_SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_ATP,
        CHEBI_GTP,
        CHEBI_UTP,
        CHEBI_CTP,
    }
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_H_PLUS,
        CHEBI_ADP,
        CHEBI_GDP,
        CHEBI_UDP,
        CHEBI_CDP,
        CHEBI_UMP,
        CHEBI_PHOSPHATE,
        CHEBI_DIPHOSPHATE,
        CHEBI_POLYPHOSPHATE,
    }
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_O2,
        CHEBI_H2O2,
    }
    CYCLIC_AMIDE_PATTERN = Chem.MolFromSmarts("[N;R][C;R](=[O])[C,N;R]")
    CYCLIC_AMIDINE_PATTERN = Chem.MolFromSmarts("[N;R][C;R](=[N;R])[C,N;R]")
    AMIDINE_PATTERN = Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3,NX2]")
    LINEAR_AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;!R]")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,X1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for water-dependent hydrolysis of cyclic amides."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not cyclic amide hydrolysis",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Cyclic amide hydrolases require water as a reactant",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Redox or oxygen chemistry indicates a different reaction class",
            )

        left_reactive = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        left_cyclic_amides = self._pattern_count(left_reactive, self.CYCLIC_AMIDE_PATTERN)
        right_cyclic_amides = self._pattern_count(right_reactive, self.CYCLIC_AMIDE_PATTERN)
        if left_cyclic_amides <= right_cyclic_amides:
            return ClassificationResult(
                is_member=False,
                explanation="Cyclic amide count does not decrease across the reaction",
            )

        if (
            self._pattern_count(left_reactive, self.CYCLIC_AMIDINE_PATTERN) > 0
            or self._pattern_count(left_reactive, self.AMIDINE_PATTERN) > 0
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Cyclic amidine substrates belong to a different hydrolase branch",
            )

        left_carboxyls = self._pattern_count(left_reactive, self.CARBOXYL_PATTERN)
        right_carboxyls = self._pattern_count(right_reactive, self.CARBOXYL_PATTERN)
        left_linear_amides = self._pattern_count(left_reactive, self.LINEAR_AMIDE_PATTERN)
        right_linear_amides = self._pattern_count(right_reactive, self.LINEAR_AMIDE_PATTERN)
        left_rings = self._ring_count(left_reactive)
        right_rings = self._ring_count(right_reactive)

        if right_carboxyls > left_carboxyls or right_linear_amides >= left_linear_amides:
            detail = "ATP-coupled ring opening" if self._has_triphosphate_coupling(reaction) else "ring opening"
            return ClassificationResult(
                is_member=True,
                explanation=f"Cyclic amide hydrolase: water-dependent {detail} of a lactam-like substrate",
            )

        if right_rings < left_rings:
            return ClassificationResult(
                is_member=True,
                explanation="Cyclic amide hydrolase: ring count decreases after hydrolysis",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No cyclic-amide hydrolysis signature detected",
        )

    @classmethod
    def _has_triphosphate_coupling(cls, reaction: Reaction) -> bool:
        """Check for ATP-like activation coupled to hydrolysis."""
        has_ntp = any(
            participant.chebi_id in {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}
            for participant in reaction.left_participants
        )
        has_spent = any(
            participant.chebi_id in {CHEBI_ADP, CHEBI_GDP, CHEBI_UDP, CHEBI_CDP, CHEBI_UMP}
            for participant in reaction.right_participants
        )
        has_phosphate = any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, CHEBI_POLYPHOSPHATE}
            for participant in reaction.right_participants
        )
        return has_ntp and has_spent and has_phosphate

    @classmethod
    def _ring_count(cls, participants: list[Participant]) -> int:
        """Count rings across substantive participants."""
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += mol.GetRingInfo().NumRings() * max(participant.count, 1)
        return total

    @classmethod
    def _pattern_count(cls, participants: list[Participant], pattern: Chem.Mol | None) -> int:
        """Count SMARTS matches across a participant list."""
        if pattern is None:
            return 0
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total
