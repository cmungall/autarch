"""Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a linear amide."""

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


class HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInLinearAmides(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds in linear amides

    Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a linear amide.
    """

    GO_ID = "GO:0016811"
    EC_NUMBER_PREFIX = "3.5.1.-"
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
    LINEAR_AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;!R]")
    PEPTIDE_PATTERN = Chem.MolFromSmarts("[NX3][CH]([#6])[CX3](=[OX1])[NX3][CH]([#6])")
    FREE_AMINE_PATTERN = Chem.MolFromSmarts("[N;X3,X4+;H1,H2,H3+;!$(N[C,S,P]=O)]")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,X1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of linear amides."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not linear amide hydrolysis",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Linear amide hydrolases require water as a reactant",
            )

        if self._has_excluded_cofactor(reaction):
            return ClassificationResult(
                is_member=False,
                explanation="Cofactor-driven chemistry indicates a different reaction class",
            )

        left_reactive = [
            participant for participant in reaction.left_participants if participant.chebi_id != CHEBI_H2O
        ]
        right_reactive = reaction.right_participants
        if any(participant.get_mol() is None for participant in left_reactive + right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for substantive participants",
            )

        left_amides = self._pattern_count(left_reactive, self.LINEAR_AMIDE_PATTERN)
        right_amides = self._pattern_count(right_reactive, self.LINEAR_AMIDE_PATTERN)
        if left_amides <= right_amides:
            return ClassificationResult(
                is_member=False,
                explanation="Linear amide count does not decrease across the reaction",
            )

        has_ammonia_product = self._has_ammonia_product(right_reactive)
        if (
            self._pattern_count(left_reactive, self.PEPTIDE_PATTERN)
            > self._pattern_count(right_reactive, self.PEPTIDE_PATTERN)
            and not has_ammonia_product
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Peptide-bond hydrolysis is excluded from the linear-amide subclass",
            )

        if has_ammonia_product:
            return ClassificationResult(
                is_member=True,
                explanation="Linear amide hydrolase: water-dependent amide cleavage with ammonia release",
            )

        left_amines = self._pattern_count(left_reactive, self.FREE_AMINE_PATTERN)
        right_amines = self._pattern_count(right_reactive, self.FREE_AMINE_PATTERN)
        left_carboxyls = self._pattern_count(left_reactive, self.CARBOXYL_PATTERN)
        right_carboxyls = self._pattern_count(right_reactive, self.CARBOXYL_PATTERN)
        if right_amines > left_amines and right_carboxyls > left_carboxyls:
            return ClassificationResult(
                is_member=True,
                explanation="Linear amide hydrolase: deacylation exposes a free amine and carboxylate product",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No linear-amide hydrolysis signature detected",
        )

    @classmethod
    def _has_ammonia_product(cls, participants: list[Participant]) -> bool:
        """Detect ammonia or ammonium release."""
        for participant in participants:
            if participant.chebi_id in {CHEBI_NH3, CHEBI_NH4}:
                return True
            mol = participant.get_mol()
            if mol is None:
                continue
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            if carbon_count == 0 and nitrogen_count == 1:
                return True
        return False

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

    @classmethod
    def _has_excluded_cofactor(cls, reaction: Reaction) -> bool:
        """Reject obvious non-hydrolase cofactors."""
        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        return bool(chebis & cls.EXCLUDED_CHEBIS)
