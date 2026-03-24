"""phosphotransferase activity, alcohol group as acceptor.

Catalysis of the transfer of a phosphorus-containing group from one compound
(donor) to an alcohol group (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UTP,
)
from autarch.ontology.transferase import Transferase


class PhosphotransferaseAlcoholGroupAsAcceptor(Transferase):
    """phosphotransferase activity, alcohol group as acceptor.

    Catalysis of the transfer of a phosphorus-containing group from one
    compound (donor) to an alcohol group (acceptor).
    """

    GO_ID = "GO:0016773"
    EC_NUMBER_PREFIX = "2.7.1.-"
    COFACTOR_CHEBIS = {
        CHEBI_ATP,
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_GTP,
        CHEBI_GDP,
        CHEBI_CTP,
        CHEBI_CDP,
        CHEBI_UTP,
        CHEBI_UDP,
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }
    DONOR_PRODUCT_PAIRS = (
        (CHEBI_ATP, CHEBI_ADP),
        (CHEBI_ADP, CHEBI_AMP),
        (CHEBI_GTP, CHEBI_GDP),
        (CHEBI_CTP, CHEBI_CDP),
        (CHEBI_UTP, CHEBI_UDP),
    )
    PHOSPHORIC_ANHYDRIDE_PATTERN = Chem.MolFromSmarts("[P][O][P]")
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])O[P]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for phosphoryl transfer onto an alcohol-bearing acceptor."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not alcohol-group phosphotransferase chemistry",
            )

        has_donor_pair = any(
            any(participant.chebi_id == donor for participant in reaction.left_participants)
            and any(participant.chebi_id == product for participant in reaction.right_participants)
            for donor, product in self.DONOR_PRODUCT_PAIRS
        )
        if not has_donor_pair:
            return ClassificationResult(
                is_member=False,
                explanation="No nucleotide phosphate-donor -> nucleotide-product coupling pattern",
            )

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.COFACTOR_CHEBIS
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.COFACTOR_CHEBIS
        ]
        if not left_core or not right_core:
            return ClassificationResult(
                is_member=False,
                explanation="No non-cofactor substrates/products to evaluate",
            )

        has_alcohol_acceptor = any(self._is_alcohol_acceptor(participant) for participant in left_core)
        if not has_alcohol_acceptor:
            return ClassificationResult(
                is_member=False,
                explanation="No alcohol-bearing acceptor detected on the substrate side",
            )

        if any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, "CHEBI:16838", "CHEBI:33019"}
            for participant in reaction.right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Free phosphate release indicates activation chemistry rather than direct alcohol-group phosphorylation",
            )

        left_anhydrides = sum(self._phosphoric_anhydride_count(participant) for participant in left_core)
        right_anhydrides = sum(self._phosphoric_anhydride_count(participant) for participant in right_core)
        if right_anhydrides > left_anhydrides:
            return ClassificationResult(
                is_member=False,
                explanation="Formation of a new phosphoric anhydride indicates phosphate-group acceptor chemistry, not alcohol-group phosphorylation",
            )

        if any(self._is_acyl_phosphate(participant) for participant in right_core):
            return ClassificationResult(
                is_member=False,
                explanation="Acyl-phosphate formation indicates carboxyl-group acceptor chemistry, not alcohol-group phosphorylation",
            )

        left_total_phosphorus = sum(self._phosphorus_count(participant) for participant in left_core)
        right_total_phosphorus = sum(self._phosphorus_count(participant) for participant in right_core)
        has_phosphorylated_product = any(
            self._phosphorus_count(participant) > 0 for participant in right_core
        )
        if not has_phosphorylated_product or right_total_phosphorus <= left_total_phosphorus:
            return ClassificationResult(
                is_member=False,
                explanation="No newly phosphorylated product detected outside the nucleotide cofactors",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Alcohol-group phosphotransferase: nucleotide-dependent phosphorylation of an alcohol acceptor",
        )

    @staticmethod
    def _is_alcohol_acceptor(participant: Participant) -> bool:
        return bool(
            participant.smiles
            and participant.has_moiety(Moiety.HYDROXYL)
            and not participant.is_thioester()
        )

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        """Count phosphorus atoms in a participant."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _phosphoric_anhydride_count(cls, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None or cls.PHOSPHORIC_ANHYDRIDE_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(cls.PHOSPHORIC_ANHYDRIDE_PATTERN))

    @classmethod
    def _is_acyl_phosphate(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None or cls.ACYL_PHOSPHATE_PATTERN is None:
            return False
        return mol.HasSubstructMatch(cls.ACYL_PHOSPHATE_PATTERN)
