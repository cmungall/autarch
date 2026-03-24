"""hydrolase acting on carbon-nitrogen but not peptide bonds in cyclic amidines.

Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
cyclic amidine, a compound of the form R-C(=NH)-NH2.
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


class HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmidines(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds in cyclic amidines.

    Catalysis of the hydrolysis of any non-peptide carbon-nitrogen bond in a
    cyclic amidine, a compound of the form R-C(=NH)-NH2.
    """

    GO_ID = "GO:0016814"
    EC_NUMBER_PREFIX = "3.5.4.-"
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
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")
    LINEAR_AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;!R]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolytic deamination within the EC 3.5.4 branch."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not cyclic amidine hydrolysis",
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
                explanation="Cyclic amidine hydrolases require water as a reactant",
            )

        if not any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="Cyclic amidine hydrolases release ammonia or ammonium",
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
                explanation="Cyclic amidine hydrolysis requires one substantive substrate and one substantive product",
            )

        substrate = left_reactive[0]
        product = right_reactive[0]
        substrate_mol = substrate.get_mol()
        product_mol = product.get_mol()
        if substrate_mol is None or product_mol is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for cyclic amidine substrates or products",
            )

        if substrate_mol.GetRingInfo().NumRings() == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Substrate lacks the cyclic scaffold expected for EC 3.5.4 chemistry",
            )
        if self._ring_nitrogen_count(substrate) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Substrate lacks ring nitrogen atoms characteristic of cyclic amidine deaminase chemistry",
            )
        if self._pattern_count(substrate, self.LINEAR_AMIDE_PATTERN) > 0:
            return ClassificationResult(
                is_member=False,
                explanation="Linear amide hydrolysis belongs to a different hydrolase branch",
            )

        if self._carbon_count(product) != self._carbon_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolytic deamination should preserve the carbon skeleton",
            )

        if self._nitrogen_count(product) != self._nitrogen_count(substrate) - 1:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not lose exactly one nitrogen atom",
            )

        if self._oxygen_count(product) != self._oxygen_count(substrate) + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not gain exactly one oxygen atom from hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Cyclic amidine hydrolase: water-driven deamination of a cyclic nitrogen scaffold",
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
    def _ring_nitrogen_count(participant: Participant) -> int:
        """Count nitrogens that are part of a ring system."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(
            1
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 7 and atom.IsInRing()
        )

    @staticmethod
    def _pattern_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count substructure matches for a participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))
