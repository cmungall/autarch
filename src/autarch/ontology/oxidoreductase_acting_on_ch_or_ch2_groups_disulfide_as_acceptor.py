"""oxidoreductase acting on CH or CH2 groups, disulfide as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH2 group acts
as a hydrogen or electron donor and reduces a disulfide group.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_FMN,
    CHEBI_FMNH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"
PROTEIN_DISULFIDE_CHEBIS = {"CHEBI:50058"}
PROTEIN_DITHIOL_CHEBIS = {"CHEBI:29950"}
PHYLLOQUINOL_LIKE_CHEBIS = {"CHEBI:28433", "CHEBI:18298"}
PHYLLOQUINONE_LIKE_CHEBIS = {"CHEBI:18067", "CHEBI:15759"}


class OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on CH or CH2 groups, disulfide as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH2 group
    acts as a hydrogen or electron donor and reduces a disulfide group.
    """

    GO_ID = "GO:0016728"
    EC_NUMBER_PREFIX = "1.17.4.-"
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_FERREDOXIN_OXIDIZED,
        CHEBI_FERREDOXIN_REDUCED,
    }
    SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    DISULFIDE_PATTERN = Chem.MolFromSmarts("[S;X2][S;X2]")
    THIOL_PATTERN = Chem.MolFromSmarts("[S;X2H1,X1-;!$(S(=O));!$(SS)]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")
    SULFUR_OXYGEN_PATTERN = Chem.MolFromSmarts("[S;X3,X4](=[O;X1])")
    EPOXIDE_PATTERN = Chem.MolFromSmarts("[OX2r3]1[#6r3][#6r3]1")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for CH/CH2 redox coupled to disulfide reduction."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH/CH2 disulfide redox",
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
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide, flavin, ferredoxin, peroxide, or dioxygen cofactors indicate another CH/CH2 oxidoreductase branch",
            )

        if not self._has_disulfide_acceptor_signature(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No disulfide acceptor is reduced in this reaction orientation",
            )

        if self._has_sulfur_donor_signature(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Partner sulfur oxidation indicates sulfur-group donor chemistry rather than CH/CH2 donor chemistry",
            )

        left_carbon = self._organic_non_sulfur(left_participants)
        right_carbon = self._organic_non_sulfur(right_participants)
        if not left_carbon or not right_carbon:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive non-sulfur carbon donor chemistry remains after removing sulfur acceptor participants",
            )

        if all(
            any(left_participant.is_same_molecule(right_participant) for right_participant in right_carbon)
            for left_participant in left_carbon
        ):
            return ClassificationResult(
                is_member=False,
                explanation="No change in non-sulfur donor scaffold detected",
            )

        if self._has_deoxyribonucleotide_oxidation(left_carbon, right_carbon):
            return ClassificationResult(
                is_member=True,
                explanation="CH/CH2 oxidoreductase: a disulfide acceptor is reduced as a deoxyribonucleotide is oxidized to a ribonucleotide",
            )

        if self._has_phylloquinol_oxidation(left_participants, right_participants):
            return ClassificationResult(
                is_member=True,
                explanation="CH/CH2 oxidoreductase: a disulfide acceptor is reduced as a phylloquinol-like donor is oxidized",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Non-sulfur donor does not match the deoxyribonucleotide or phylloquinol-like CH/CH2 redox signatures in this branch",
        )

    def _has_disulfide_acceptor_signature(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        left_chebis = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if left_chebis & PROTEIN_DISULFIDE_CHEBIS and right_chebis & PROTEIN_DITHIOL_CHEBIS:
            return True
        return self._count_pattern(left_participants, self.DISULFIDE_PATTERN) >= 1 and self._count_pattern(
            right_participants, self.THIOL_PATTERN
        ) >= 2

    def _has_sulfur_donor_signature(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        left_thiols = self._count_pattern(left_participants, self.THIOL_PATTERN)
        right_disulfides = self._count_pattern(right_participants, self.DISULFIDE_PATTERN)
        if left_thiols >= 2 and right_disulfides >= 1:
            return True

        return self._count_pattern(right_participants, self.SULFUR_OXYGEN_PATTERN) > self._count_pattern(
            left_participants, self.SULFUR_OXYGEN_PATTERN
        )

    def _organic_non_sulfur(self, participants: list[Participant]) -> list[Participant]:
        organic: list[Participant] = []
        for participant in participants:
            if participant.chebi_id in self.SPECTATOR_CHEBIS:
                continue
            mol = participant.get_mol()
            if mol is None:
                continue
            atom_numbers = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
            if 6 not in atom_numbers or 16 in atom_numbers:
                continue
            organic.append(participant)
        return organic

    def _has_deoxyribonucleotide_oxidation(
        self,
        left_carbon: list[Participant],
        right_carbon: list[Participant],
    ) -> bool:
        for left_participant in left_carbon:
            if self._phosphorus_count(left_participant) < 2:
                continue
            for right_participant in right_carbon:
                if self._phosphorus_count(left_participant) != self._phosphorus_count(right_participant):
                    continue
                if self._carbon_count(left_participant) != self._carbon_count(right_participant):
                    continue
                if self._oxygen_count([right_participant]) != self._oxygen_count([left_participant]) + 1:
                    continue
                return True
        return False

    def _has_phylloquinol_oxidation(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        left_chebis = {participant.chebi_id for participant in left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in right_participants if participant.chebi_id}
        if not (left_chebis & PHYLLOQUINOL_LIKE_CHEBIS and right_chebis & PHYLLOQUINONE_LIKE_CHEBIS):
            return False
        return self._count_pattern(right_participants, self.EPOXIDE_PATTERN) >= self._count_pattern(
            left_participants, self.EPOXIDE_PATTERN
        )

    @staticmethod
    def _count_pattern(participants: list[Participant], pattern: Chem.Mol | None) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None or pattern is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total

    @staticmethod
    def _oxygen_count(participants: list[Participant]) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8) * max(
                participant.count, 1
            )
        return total

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
