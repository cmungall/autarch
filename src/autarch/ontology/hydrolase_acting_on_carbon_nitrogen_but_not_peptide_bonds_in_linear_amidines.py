"""hydrolase acting on carbon-nitrogen but not peptide bonds in linear amidines.

Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
linear amidine, a compound of the form R-C(=NH)-NH2.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.hydrolase import Hydrolase


class HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmidines(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds in linear
    amidines.

    Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
    linear amidine, a compound of the form R-C(=NH)-NH2.
    """

    GO_ID = "GO:0016813"
    EC_NUMBER_PREFIX = "3.5.3.-"
    EXCLUDED_CHEBIS = {
        CHEBI_ATP,
        CHEBI_COA,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_O2,
    }
    LINEAR_AMIDINE_PATTERN = Chem.MolFromSmarts("[CX3;!R](=[NX2,NX3+])[NX3;!R]")
    LINEAR_AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;!R]")
    UREA_PATTERN = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of a linear amidine."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not linear amidine hydrolysis",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Linear amidine hydrolases require water as a reactant",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Cofactor-driven chemistry indicates a different reaction class",
            )

        left_reactive = [
            participant for participant in reaction.left_participants if participant.chebi_id != CHEBI_H2O
        ]
        right_reactive = reaction.right_participants
        if any(
            participant.get_mol() is None
            for participant in left_reactive + [p for p in right_reactive if p.chebi_id not in {CHEBI_NH3, CHEBI_NH4}]
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for substantive participants",
            )

        left_amidines = self._pattern_count(left_reactive, self.LINEAR_AMIDINE_PATTERN)
        right_amidines = self._pattern_count(right_reactive, self.LINEAR_AMIDINE_PATTERN)
        if left_amidines <= right_amidines:
            return ClassificationResult(
                is_member=False,
                explanation="Linear amidine count does not decrease across the reaction",
            )

        if self._pattern_count(right_reactive, self.UREA_PATTERN) > 0:
            return ClassificationResult(
                is_member=True,
                explanation="Linear amidine hydrolase: amidine hydrolysis yields a urea-containing product",
            )

        left_amides = self._pattern_count(left_reactive, self.LINEAR_AMIDE_PATTERN)
        right_amides = self._pattern_count(right_reactive, self.LINEAR_AMIDE_PATTERN)
        if right_amides > left_amides:
            return ClassificationResult(
                is_member=True,
                explanation="Linear amidine hydrolase: amidine hydrolysis yields an amide-containing product",
            )

        has_ammonia = any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right_reactive)
        if has_ammonia and right_amides >= left_amides:
            return ClassificationResult(
                is_member=True,
                explanation="Linear amidine hydrolase: amidine hydrolysis with ammonia release",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No linear amidine hydrolysis signature detected",
        )

    @staticmethod
    def _pattern_count(participants: list[Participant], pattern: Chem.Mol | None) -> int:
        if pattern is None:
            return 0
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total
