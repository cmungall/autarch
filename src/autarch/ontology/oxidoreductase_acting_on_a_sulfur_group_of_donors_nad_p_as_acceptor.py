"""oxidoreductase activity, acting on a sulfur group of donors, NAD(P) as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a sulfur-containing
group acts as a hydrogen or electron donor and reduces NAD or NADP.
"""

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase_acting_on_a_sulfur_group_of_donors import (
    OxidoreductaseActingOnASulfurGroupOfDonors,
)


class OxidoreductaseActingOnASulfurGroupOfDonorsNADPAsAcceptor(
    OxidoreductaseActingOnASulfurGroupOfDonors
):
    """oxidoreductase activity, acting on a sulfur group of donors, NAD(P) as
    acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a
    sulfur-containing group acts as a hydrogen or electron donor and reduces NAD
    or NADP.
    """

    GO_ID = "GO:0016668"
    EC_NUMBER_PREFIX = "1.8.1.-"
    OXIDIZED_NICOTINAMIDES = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    REDUCED_NICOTINAMIDES = {CHEBI_NADH, CHEBI_NADPH}
    LEFT_SPECTATOR_CHEBIS = OXIDIZED_NICOTINAMIDES | {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = REDUCED_NICOTINAMIDES | {CHEBI_H2O, CHEBI_H_PLUS}
    EXCLUDED_EXTERNAL_CHEBIS = {CHEBI_O2, CHEBI_H2O2, CHEBI_FAD, CHEBI_FADH2}
    THIOL_PATTERN = Chem.MolFromSmarts("[S;X2H1,X1-;!$(S(=O));!$(SS)]")
    DISULFIDE_PATTERN = Chem.MolFromSmarts("[S;X2][S;X2]")
    SULFUR_OXYGEN_PATTERN = Chem.MolFromSmarts("[S;X3,X4](=[O;X1])")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for nicotinamide-linked sulfur-group redox chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sulfur-group oxidoreductase chemistry",
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
                explanation="Sulfur-group oxidoreductase detected in the reverse nicotinamide direction",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No sulfur-group NAD(P)-linked redox signature detected",
        )

    def _check_direction(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        if not self._has_nicotinamide_pair(left, right):
            return ClassificationResult(
                is_member=False,
                explanation="Requires NAD(+) or NADP(+) reduction across the reaction",
            )

        chebis = {
            participant.chebi_id
            for participant in left + right
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_EXTERNAL_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External oxygen or flavin cofactors indicate a different oxidoreductase branch",
            )

        left_reactive = self._reactive_participants(left, self.LEFT_SPECTATOR_CHEBIS)
        right_reactive = self._reactive_participants(right, self.RIGHT_SPECTATOR_CHEBIS)
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive sulfur-containing participants remain after removing nicotinamides and solvent spectators",
            )

        if self._has_thiol_disulfide_redox(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur-group oxidoreductase: nicotinamide-linked thiol/disulfide redox",
            )

        if self._has_sulfur_oxygenation_redox(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur-group oxidoreductase: nicotinamide-linked sulfur oxyanion redox",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No sulfur redox transformation matched the nicotinamide-coupled branch",
        )

    def _has_nicotinamide_pair(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        left_oxidized = self._count_chebis(left, self.OXIDIZED_NICOTINAMIDES)
        right_reduced = self._count_chebis(right, self.REDUCED_NICOTINAMIDES)
        left_reduced = self._count_chebis(left, self.REDUCED_NICOTINAMIDES)
        right_oxidized = self._count_chebis(right, self.OXIDIZED_NICOTINAMIDES)
        return (
            left_oxidized >= 1
            and left_oxidized == right_reduced
            and left_reduced == 0
            and right_oxidized == 0
        )

    @staticmethod
    def _count_chebis(participants: list[Participant], chebis: set[str]) -> int:
        return sum(
            max(participant.count, 1)
            for participant in participants
            if participant.chebi_id in chebis
        )

    @staticmethod
    def _reactive_participants(
        participants: list[Participant],
        spectator_chebis: set[str],
    ) -> list[Participant]:
        return [
            participant
            for participant in participants
            if participant.chebi_id not in spectator_chebis and participant.get_mol() is not None
        ]

    @classmethod
    def _has_thiol_disulfide_redox(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        left_thiols = sum(
            cls._match_count(participant, cls.THIOL_PATTERN) * max(participant.count, 1)
            for participant in left_reactive
        )
        right_disulfides = sum(
            cls._match_count(participant, cls.DISULFIDE_PATTERN) * max(participant.count, 1)
            for participant in right_reactive
        )
        left_sulfur = sum(
            cls._sulfur_count(participant) * max(participant.count, 1)
            for participant in left_reactive
        )
        right_sulfur = sum(
            cls._sulfur_count(participant) * max(participant.count, 1)
            for participant in right_reactive
        )
        return left_thiols >= 2 and right_disulfides >= 1 and left_sulfur == right_sulfur

    @classmethod
    def _has_sulfur_oxygenation_redox(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        for substrate in left_reactive:
            if cls._sulfur_count(substrate) == 0:
                continue
            substrate_non_s = cls._non_sulfur_heavy_atom_counts(substrate)
            substrate_c = cls._carbon_count(substrate)
            substrate_so = cls._sulfur_oxygen_count(substrate)
            for product in right_reactive:
                if cls._sulfur_count(product) != cls._sulfur_count(substrate):
                    continue
                if substrate_non_s != cls._non_sulfur_heavy_atom_counts(product):
                    if not (substrate_c == 0 and cls._carbon_count(product) == 0):
                        continue
                if cls._sulfur_oxygen_count(product) > substrate_so:
                    return True
        return False

    @staticmethod
    def _match_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))

    @staticmethod
    def _sulfur_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @classmethod
    def _non_sulfur_heavy_atom_counts(cls, participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() not in {1, 16}
        )

    @classmethod
    def _sulfur_oxygen_count(cls, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 16:
                continue
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() >= 1.0:
                    count += 1
        count += cls._match_count(participant, cls.SULFUR_OXYGEN_PATTERN)
        return count
