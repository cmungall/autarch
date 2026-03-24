"""sulfuric ester hydrolase.

Catalysis of the reaction: RSO-R' + H2O = RSOOH + R'H. This reaction is the
hydrolysis of a sulfuric ester bond, an ester formed from sulfuric acid,
O=SO(OH)2.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
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


class SulfuricEsterHydrolase(Hydrolase):
    """sulfuric ester hydrolase.

    Catalysis of the reaction: RSO-R' + H2O = RSOOH + R'H. This reaction is
    the hydrolysis of a sulfuric ester bond, an ester formed from sulfuric acid,
    O=SO(OH)2.
    """

    GO_ID = "GO:0008484"
    EC_NUMBER_PREFIX = "3.1.6.-"
    EXCLUDED_CHEBIS = {
        CHEBI_ATP,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_O2,
    }
    SULFATE_ESTER_PATTERN = Chem.MolFromSmarts("[#6][O][S](=O)(=O)[O-,$([OH])]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for hydrolysis of an O-sulfate ester."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sulfuric ester hydrolysis",
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
        """Evaluate one reaction orientation."""
        if not any(participant.chebi_id == CHEBI_H2O for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="Sulfuric ester hydrolases require water as a reactant",
            )

        chebis = {
            participant.chebi_id
            for participant in left + right
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide, oxygen, or redox cofactors indicate different chemistry",
            )

        left_reactive = [
            participant
            for participant in left
            if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_reactive = [
            participant
            for participant in right
            if participant.chebi_id != CHEBI_H_PLUS
        ]
        if any(participant.get_mol() is None for participant in left_reactive + right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for sulfuric ester hydrolysis",
            )

        left_sulfate_esters = [participant for participant in left_reactive if self._is_sulfate_ester(participant)]
        right_sulfate_esters = [participant for participant in right_reactive if self._is_sulfate_ester(participant)]
        if not left_sulfate_esters or len(right_sulfate_esters) >= len(left_sulfate_esters):
            return ClassificationResult(
                is_member=False,
                explanation="No sulfuric ester bond is hydrolyzed",
            )

        if not any(self._is_inorganic_sulfate(participant) for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No inorganic sulfate product detected",
            )

        if any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia-releasing sulfate elimination indicates different chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sulfuric ester hydrolase: water-dependent cleavage of an O-sulfate ester to sulfate and a desulfated product",
        )

    @classmethod
    def _is_sulfate_ester(cls, participant: Participant) -> bool:
        if participant.is_sulfamate():
            return False
        mol = participant.get_mol()
        return mol is not None and cls.SULFATE_ESTER_PATTERN is not None and mol.HasSubstructMatch(
            cls.SULFATE_ESTER_PATTERN
        )

    @staticmethod
    def _is_inorganic_sulfate(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        return carbon_count == 0 and sulfur_count == 1 and oxygen_count >= 4
