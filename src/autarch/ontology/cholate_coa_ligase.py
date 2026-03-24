"""cholate-CoA ligase activity.

Catalysis of ATP-dependent formation of a bile acid CoA thioester from cholate
or a closely related polyhydroxylated steroid acid.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.ligase_forming_carbon_sulfur_bonds import LigaseFormingCarbonSulfurBonds


class CholateCoALigase(LigaseFormingCarbonSulfurBonds):
    """cholate-CoA ligase activity.

    Catalysis of ATP-dependent formation of a bile acid CoA thioester from
    cholate or a closely related polyhydroxylated steroid acid.
    """

    GO_ID = "GO:0047747"
    EC_NUMBER_PREFIX = "6.2.1.7"
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[OX2H]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP-dependent CoA ligation of a bile acid-like substrate."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_acids = [
            participant
            for participant in reaction.left_participants
            if self._is_bile_acid_like(participant, require_thioester=False)
        ]
        right_products = [
            participant
            for participant in reaction.right_participants
            if self._is_bile_acid_like(participant, require_thioester=True)
        ]
        if not left_acids:
            return ClassificationResult(
                is_member=False,
                explanation="No bile acid-like acid substrate detected",
            )
        if not right_products:
            return ClassificationResult(
                is_member=False,
                explanation="No bile acid-like CoA thioester product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Cholate-CoA ligase: ATP-dependent CoA ligation of a polyhydroxylated steroid acid",
        )

    @classmethod
    def _is_bile_acid_like(cls, participant: Participant, *, require_thioester: bool) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        if require_thioester and not participant.is_thioester():
            return False
        if not require_thioester and (participant.is_thioester() or not participant.is_carboxylic_acid()):
            return False
        ring_count = mol.GetRingInfo().NumRings()
        hydroxyl_count = len(mol.GetSubstructMatches(cls.HYDROXYL_PATTERN)) if cls.HYDROXYL_PATTERN else 0
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        return ring_count >= 4 and hydroxyl_count >= 3 and carbon_count >= 24
