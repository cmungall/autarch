"""intramolecular oxidoreductase interconverting aldoses and ketoses.

Catalysis of an oxidation-reduction (redox) reaction in which the hydrogen donor
and acceptor, which is an aldose or a ketose, are the same molecule, and no
oxidized product appears.
"""

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.isomerase import Isomerase


class IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses(Isomerase):
    """intramolecular oxidoreductase interconverting aldoses and ketoses.

    Catalysis of an oxidation-reduction (redox) reaction in which the hydrogen
    donor and acceptor, which is an aldose or a ketose, are the same molecule,
    and no oxidized product appears.
    """

    GO_ID = "GO:0016861"
    EC_NUMBER_PREFIX = "5.3.1.-"
    EXCLUDED_CHEBIS = {
        CHEBI_H2O,
        CHEBI_O2,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_ATP,
        CHEBI_ADP,
        CHEBI_CO2,
    }
    ALDEHYDE_PATTERN = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")
    KETONE_PATTERN = Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]")
    TERMINAL_HYDROXYMETHYL_PATTERN = Chem.MolFromSmarts("[CH2X4][OX2,PX4,SX4]")
    RING_OXYGEN_PATTERN = Chem.MolFromSmarts("[O;R]")
    HYDROXYLATED_CARBON_PATTERN = Chem.MolFromSmarts("[CX4][OX2H,OX1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for aldose-ketose interconversion without external cofactors."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not intramolecular aldose-ketose isomerization",
            )

        if len(reaction.left_participants) != 1 or len(reaction.right_participants) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Aldose-ketose isomerases operate on a single substrate-product pair",
            )

        if self._has_excluded_cofactor(reaction):
            return ClassificationResult(
                is_member=False,
                explanation="External cofactors indicate a different reaction class",
            )

        substrate = reaction.left_participants[0]
        product = reaction.right_participants[0]
        if substrate.get_mol() is None or product.get_mol() is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for substrate or product",
            )

        if not self._same_heavy_atom_composition(substrate, product):
            return ClassificationResult(
                is_member=False,
                explanation="Aldose-ketose isomerization conserves heavy-atom composition",
            )

        if self._is_aldose_ketose_pair(substrate, product):
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: aldose-ketose interconversion without external cofactors",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No aldose-ketose structural interconversion detected",
        )

    @classmethod
    def _is_aldose_ketose_pair(cls, left: Participant, right: Participant) -> bool:
        """Accept either aldose -> ketose or ketose -> aldose orientation."""
        return (
            cls._is_aldose_like(left) and cls._is_ketose_like(right)
        ) or (
            cls._is_ketose_like(left) and cls._is_aldose_like(right)
        )

    @classmethod
    def _is_aldose_like(cls, participant: Participant) -> bool:
        """Detect aldose-like chemistry, including cyclic hemiacetal forms."""
        if cls._has_pattern(participant, cls.ALDEHYDE_PATTERN):
            return True
        ring_oxygen_count = cls._match_count(participant, cls.RING_OXYGEN_PATTERN)
        terminal_count = cls._match_count(participant, cls.TERMINAL_HYDROXYMETHYL_PATTERN)
        hydroxylated_carbons = cls._match_count(participant, cls.HYDROXYLATED_CARBON_PATTERN)
        return (
            cls._oxygen_count(participant) >= 3
            and ring_oxygen_count >= 1
            and hydroxylated_carbons >= 3
            and not (ring_oxygen_count == 1 and terminal_count >= 2)
            and not cls._has_pattern(participant, cls.KETONE_PATTERN)
        )

    @classmethod
    def _is_ketose_like(cls, participant: Participant) -> bool:
        """Detect ketose-like chemistry, including cyclic ketose forms."""
        if cls._has_pattern(participant, cls.KETONE_PATTERN):
            return True
        return (
            cls._oxygen_count(participant) >= 3
            and cls._match_count(participant, cls.RING_OXYGEN_PATTERN) == 1
            and cls._match_count(participant, cls.HYDROXYLATED_CARBON_PATTERN) >= 3
            and cls._match_count(participant, cls.TERMINAL_HYDROXYMETHYL_PATTERN) >= 2
            and not cls._has_pattern(participant, cls.ALDEHYDE_PATTERN)
        )

    @classmethod
    def _same_heavy_atom_composition(cls, left: Participant, right: Participant) -> bool:
        """Compare heavy-atom composition."""
        return cls._heavy_atom_counts(left) == cls._heavy_atom_counts(right)

    @classmethod
    def _heavy_atom_counts(cls, participant: Participant) -> Counter[int]:
        """Count atoms by atomic number, excluding hydrogen."""
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() != 1
        )

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        """Count oxygen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @staticmethod
    def _has_excluded_cofactor(reaction: Reaction) -> bool:
        """Reject reactions with external cofactors or obvious non-isomerase participants."""
        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        return bool(chebis & IntramolecularOxidoreductaseInterconvertingAldosesAndKetoses.EXCLUDED_CHEBIS)

    @staticmethod
    def _match_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count SMARTS matches in one participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))

    @classmethod
    def _has_pattern(cls, participant: Participant, pattern: Chem.Mol | None) -> bool:
        """Check whether a participant matches a SMARTS pattern."""
        return cls._match_count(participant, pattern) > 0
