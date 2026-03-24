"""acid-amino acid ligase.

Catalysis of the ligation of an acid to an amino acid via a carbon-nitrogen bond,
with the concomitant hydrolysis of the diphosphate bond in ATP or a similar
triphosphate.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_COA,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.ligase import Ligase

CHEBI_TRNA = "CHEBI:78442"


class AcidAminoAcidLigase(Ligase):
    """acid-amino acid ligase.

    Catalysis of the ligation of an acid to an amino acid via a carbon-nitrogen
    bond, with the concomitant hydrolysis of the diphosphate bond in ATP or a
    similar triphosphate.
    """

    GO_ID = "GO:0016881"
    EC_NUMBER_PREFIX = "6.3.2.-"
    NTP_REACTANTS = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}
    SPENT_NUCLEOTIDES = {
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_GDP,
        CHEBI_UDP,
        CHEBI_UMP,
        CHEBI_CDP,
        CHEBI_CMP,
    }
    LEFT_SPECTATOR_CHEBIS = NTP_REACTANTS | {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = SPENT_NUCLEOTIDES | {
        CHEBI_PHOSPHATE,
        CHEBI_DIPHOSPHATE,
        "CHEBI:16838",
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    FREE_AMINE_PATTERN = Chem.MolFromSmarts("[N;X3,X4+;H1,H2,H3+;!$(N[C,S,P]=O)]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP- or NTP-dependent acylation of a monomeric amino acid."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not acid-amino acid ligation",
            )

        forward_result = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward_result.is_member:
            return forward_result

        reverse_result = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse_result.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse_result.explanation} (reverse reaction orientation)",
            )

        return forward_result

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        if not self._has_triphosphate_coupling(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No ATP-or-similar triphosphate coupling detected",
            )

        if any(self._is_free_ammonia(participant) for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Uses free ammonia/ammonium rather than amino-acid ligation chemistry",
            )

        if any(
            participant.chebi_id in {CHEBI_COA, CHEBI_TRNA, CHEBI_O2}
            for participant in left_participants + right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="CoA, tRNA, or O2 chemistry indicates a different ligase branch",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        if len(left_reactive) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Requires at least two substantive substrates",
            )

        amino_acid_substrates = [
            participant for participant in left_reactive if self._is_amino_acid_like(participant)
        ]
        if not amino_acid_substrates:
            return ClassificationResult(
                is_member=False,
                explanation="No amino-acid-like substrate detected",
            )

        acid_substrates = [
            participant for participant in left_reactive if participant.is_carboxylic_acid()
        ]
        if not acid_substrates:
            return ClassificationResult(
                is_member=False,
                explanation="No carboxylic-acid substrate detected",
            )

        if len(acid_substrates) < 2 and not any(
            self._is_free_amine_partner(participant) for participant in left_reactive
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Acid-amino acid ligases require either two acid substrates or an amino-acid acyl donor plus a free amine acceptor",
            )

        left_amide_bonds = self._amide_bond_count(left_reactive)
        right_amide_bonds = self._amide_bond_count(right_reactive)
        if right_amide_bonds <= left_amide_bonds:
            return ClassificationResult(
                is_member=False,
                explanation="No new amide bond is formed in the principal products",
            )

        if not any(self._contains_amide(participant) for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Principal products do not retain the expected amide linkage",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Acid-amino acid ligase: triphosphate-coupled C-N bond formation "
                "between a carboxylic acid donor and an amino-acid nucleophile"
            ),
        )

    @classmethod
    def _has_triphosphate_coupling(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for ATP-or-similar triphosphate consumption."""
        has_ntp = any(participant.chebi_id in cls.NTP_REACTANTS for participant in left_participants)
        has_spent_nucleotide = any(
            participant.chebi_id in cls.SPENT_NUCLEOTIDES for participant in right_participants
        )
        has_phosphate_product = any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, "CHEBI:16838"}
            for participant in right_participants
        )
        return has_ntp and has_spent_nucleotide and has_phosphate_product

    @classmethod
    def _amide_bond_count(cls, participants: list[Participant]) -> int:
        """Count amide bonds across participants."""
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None or cls.AMIDE_PATTERN is None:
                continue
            total += len(mol.GetSubstructMatches(cls.AMIDE_PATTERN)) * max(participant.count, 1)
        return total

    @classmethod
    def _contains_amide(cls, participant: Participant) -> bool:
        """Check whether the product contains an amide linkage."""
        mol = participant.get_mol()
        return bool(mol and cls.AMIDE_PATTERN and mol.HasSubstructMatch(cls.AMIDE_PATTERN))

    @classmethod
    def _is_amino_acid_like(cls, participant: Participant) -> bool:
        """Check for a free amino acid or amino-acid-like monomer."""
        mol = participant.get_mol()
        if mol is None or cls.FREE_AMINE_PATTERN is None:
            return False
        if not participant.is_carboxylic_acid():
            return False
        if not mol.HasSubstructMatch(cls.FREE_AMINE_PATTERN):
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        return carbon_count <= 12 and cls._amide_bond_count([participant]) <= 1

    @staticmethod
    def _is_free_ammonia(participant: Participant) -> bool:
        """Check for free ammonia/ammonium even when no ChEBI ID is provided."""
        if participant.chebi_id in {CHEBI_NH3, CHEBI_NH4}:
            return True
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        return carbon_count == 0 and nitrogen_count == 1

    @classmethod
    def _is_free_amine_partner(cls, participant: Participant) -> bool:
        """Check for a non-carbohydrate amine acceptor paired with an amino-acid acyl donor."""
        mol = participant.get_mol()
        if mol is None or cls.FREE_AMINE_PATTERN is None:
            return False
        if participant.is_carboxylic_acid():
            return False
        if not mol.HasSubstructMatch(cls.FREE_AMINE_PATTERN):
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
        return carbon_count < 12 and not (oxygen_count >= 4 and phosphorus_count == 0)
