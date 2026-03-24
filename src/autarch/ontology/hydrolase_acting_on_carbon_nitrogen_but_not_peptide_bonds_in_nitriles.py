"""hydrolase acting on carbon-nitrogen but not peptide bonds in nitriles.

Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
nitrile, a compound containing the cyano radical, -CN.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_COA,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.hydrolase import Hydrolase


class HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInNitriles(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds in nitriles.

    Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
    nitrile, a compound containing the cyano radical, -CN.
    """

    GO_ID = "GO:0016815"
    EC_NUMBER_PREFIX = "3.5.5.-"
    EXCLUDED_CHEBIS = {
        CHEBI_COA,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_O2,
    }
    NITRILE_PATTERN = Chem.MolFromSmarts("[#6]#[#7]")
    THIOCYANATE_PATTERN = Chem.MolFromSmarts("[#16]-[#6]#[#7]")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[O;H1,-]")
    CARBONYL_SULFIDE_PATTERN = Chem.MolFromSmarts("[#8]=[#6]=[#16]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for nitrile hydrolysis chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not nitrile hydrolysis",
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one hydrolysis orientation."""
        if not any(participant.chebi_id == CHEBI_H2O for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="Nitrile hydrolases require water as a reactant",
            )

        if not any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="Nitrile hydrolases release ammonia or ammonium",
            )

        left_reactive = [
            participant
            for participant in left
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_reactive = [
            participant
            for participant in right
            if participant.chebi_id not in {CHEBI_NH3, CHEBI_NH4, CHEBI_H_PLUS}
        ]
        if len(left_reactive) != 1 or len(right_reactive) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Nitrile hydrolysis requires one substantive substrate and one substantive product",
            )

        substrate = left_reactive[0]
        product = right_reactive[0]
        if substrate.get_mol() is None or product.get_mol() is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for nitrile substrates or products",
            )

        if self._pattern_count(substrate, self.THIOCYANATE_PATTERN) > 0:
            if self._pattern_count(product, self.CARBONYL_SULFIDE_PATTERN) > 0:
                return ClassificationResult(
                    is_member=True,
                    explanation="Nitrile hydrolase: thiocyanate hydrolysis yields carbonyl sulfide and ammonium",
                )
            return ClassificationResult(
                is_member=False,
                explanation="Thiocyanate hydrolysis requires carbonyl sulfide as the substantive product",
            )

        if self._pattern_count(substrate, self.NITRILE_PATTERN) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Substrate lacks a nitrile carbon-nitrogen triple bond",
            )

        if self._carbon_count(product) != self._carbon_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Nitrile hydrolysis should preserve the carbon skeleton",
            )

        if self._nitrogen_count(product) > self._nitrogen_count(substrate) - 1:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not lose the nitrile nitrogen atom",
            )

        if self._oxygen_count(product) < self._oxygen_count(substrate) + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not gain oxygen consistent with nitrile hydration",
            )

        if self._pattern_count(product, self.CARBOXYL_PATTERN) <= self._pattern_count(
            substrate, self.CARBOXYL_PATTERN
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Product does not gain carboxylate character expected from nitrile hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Nitrile hydrolase: hydration of a nitrile carbon-nitrogen triple bond to a carboxylate product",
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        """Count carbon atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        """Count nitrogen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        """Count oxygen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @staticmethod
    def _pattern_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count substructure matches for a participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))
