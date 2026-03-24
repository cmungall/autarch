"""fatty acid ligase activity.

Catalysis of the ligation of a fatty acid to an acceptor, coupled to the
hydrolysis of ATP.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.ontology.acid_thiol_ligase import AcidThiolLigase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE


class FattyAcidLigaseActivity(AcidThiolLigase):
    """fatty acid ligase activity.

    Catalysis of the ligation of a fatty acid to an acceptor, coupled to the
    hydrolysis of ATP.
    """

    GO_ID = "GO:0015645"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[C;X3](=[O;X1])[O;H1,X1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent

        if not any(participant.chebi_id in self.COA_CHEBIS for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Fatty-acid ligases require CoA as the thiol acceptor",
            )

        fatty_acids = [
            participant
            for participant in reaction.left_participants
            if self._is_fatty_acid_like(participant)
        ]
        if not fatty_acids:
            return ClassificationResult(
                is_member=False,
                explanation="No fatty-acid substrate detected",
            )

        thioester_products = [participant for participant in reaction.right_participants if participant.is_thioester()]
        if not thioester_products:
            return ClassificationResult(
                is_member=False,
                explanation="No acyl-CoA product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Fatty acid ligase activity: ATP-coupled ligation of a fatty acid to CoA",
        )

    @classmethod
    def _is_fatty_acid_like(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None or cls.CARBOXYL_PATTERN is None:
            return False
        atom_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        if any(number not in {1, 6, 8} for number in atom_numbers):
            return False
        if participant.is_thioester() or not participant.has_moiety(Moiety.CARBOXYL):
            return False
        if len(mol.GetSubstructMatches(cls.CARBOXYL_PATTERN)) != 1:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count < 6:
            return False
        return mol.GetRingInfo().NumRings() == 0
