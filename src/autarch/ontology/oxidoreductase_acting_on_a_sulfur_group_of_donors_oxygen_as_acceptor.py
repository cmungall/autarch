"""oxidoreductase activity, acting on a sulfur group of donors, oxygen as
acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a sulfur-
containing group acts as a hydrogen or electron donor and reduces oxygen.
"""

from __future__ import annotations

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
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseActingOnASulfurGroupOfDonorsOxygenAsAcceptor(Oxidoreductase):
    """oxidoreductase activity, acting on a sulfur group of donors, oxygen as
    acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a sulfur-
    containing group acts as a hydrogen or electron donor and reduces oxygen.
    """

    GO_ID = "GO:0016670"
    EC_NUMBER_PREFIX = "1.8.3.-"

    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        "CHEBI:58210",  # FMN
        "CHEBI:58307",  # FMNH2
        "CHEBI:33722",  # oxidized [4Fe-4S]
        "CHEBI:33723",  # reduced [4Fe-4S]
        "CHEBI:33737",  # oxidized ferredoxin
        "CHEBI:33738",  # reduced ferredoxin
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O2, CHEBI_H2O, CHEBI_H_PLUS}
    THIOL_PATTERN = Chem.MolFromSmarts("[S;X2H1,X1-;!$(S(=O));!$(SS)]")
    DISULFIDE_PATTERN = Chem.MolFromSmarts("[S;X2][S;X2]")
    SULFUR_OXYGEN_PATTERN = Chem.MolFromSmarts("[S;X3,X4](=[O;X1])")
    ALDEHYDE_PATTERN = Chem.MolFromSmarts("[CX3H1,H2](=[O;X1])")
    KETONE_PATTERN = Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for sulfur-donor oxidation with oxygen as the acceptor."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sulfur-group oxygenase redox chemistry",
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
        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing molecular oxygen acceptor on the substrate side",
            )
        if not any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing hydrogen peroxide product expected for oxygen-accepting sulfur redox",
            )

        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide, flavin, or iron-sulfur cofactors indicate another sulfur redox branch",
            )

        left_sulfur = self._sulfur_reactive(left_participants, self.LEFT_SPECTATOR_CHEBIS)
        right_sulfur = self._sulfur_reactive(right_participants, self.RIGHT_SPECTATOR_CHEBIS)
        right_other = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
            and self._sulfur_count(participant) == 0
        ]
        if not left_sulfur or not right_sulfur:
            return ClassificationResult(
                is_member=False,
                explanation="No sulfur-bearing donor/product pair remained after removing oxygen and solvent spectators",
            )

        if self._has_thiol_disulfide_oxidation(left_sulfur, right_sulfur):
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur-group oxidoreductase: oxygen-dependent thiol oxidation to a disulfide",
            )

        if self._has_sulfur_oxygenation(left_sulfur, right_sulfur):
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur-group oxidoreductase: oxygen-dependent sulfur oxygenation",
            )

        if self._has_oxidative_cs_cleavage(left_sulfur, right_sulfur, right_other):
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur-group oxidoreductase: oxygen-dependent oxidative C-S bond cleavage",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No sulfur-group oxidation pattern matched the oxygen-acceptor branch",
        )

    def _sulfur_reactive(
        self,
        participants: list[Participant],
        spectator_chebis: set[str],
    ) -> list[Participant]:
        return [
            participant
            for participant in participants
            if participant.chebi_id not in spectator_chebis and self._sulfur_count(participant) > 0
        ]

    def _has_thiol_disulfide_oxidation(
        self,
        left_sulfur: list[Participant],
        right_sulfur: list[Participant],
    ) -> bool:
        left_thiols = sum(self._substructure_count(p, self.THIOL_PATTERN) for p in left_sulfur)
        right_disulfides = sum(self._substructure_count(p, self.DISULFIDE_PATTERN) for p in right_sulfur)
        return left_thiols >= 2 and right_disulfides >= 1

    def _has_sulfur_oxygenation(
        self,
        left_sulfur: list[Participant],
        right_sulfur: list[Participant],
    ) -> bool:
        for substrate in left_sulfur:
            substrate_non_s = self._non_sulfur_heavy_atom_counts(substrate)
            substrate_so = self._substructure_count(substrate, self.SULFUR_OXYGEN_PATTERN)
            for product in right_sulfur:
                if self._sulfur_count(product) != self._sulfur_count(substrate):
                    continue
                if substrate_non_s != self._non_sulfur_heavy_atom_counts(product):
                    if not (
                        self._carbon_count(substrate) == 0 and self._carbon_count(product) == 0
                    ):
                        continue
                if self._substructure_count(product, self.SULFUR_OXYGEN_PATTERN) > substrate_so:
                    return True
        return False

    def _has_oxidative_cs_cleavage(
        self,
        left_sulfur: list[Participant],
        right_sulfur: list[Participant],
        right_other: list[Participant],
    ) -> bool:
        has_carbonyl_product = any(self._is_carbonyl(participant) for participant in right_other)
        if not has_carbonyl_product:
            return False
        for substrate in left_sulfur:
            substrate_s = self._sulfur_count(substrate)
            substrate_c = self._carbon_count(substrate)
            if substrate_s == 0 or substrate_c == 0:
                continue
            for product in right_sulfur:
                if self._sulfur_count(product) != substrate_s:
                    continue
                if self._carbon_count(product) < substrate_c:
                    return True
        return False

    @staticmethod
    def _substructure_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)

    @staticmethod
    def _sulfur_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16) * max(participant.count, 1)

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) * max(participant.count, 1)

    @staticmethod
    def _non_sulfur_heavy_atom_counts(participant: Participant) -> dict[int, int]:
        mol = participant.get_mol()
        if mol is None:
            return {}
        counts: dict[int, int] = {}
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num in {1, 16}:
                continue
            counts[atomic_num] = counts.get(atomic_num, 0) + 1
        return counts

    def _is_carbonyl(self, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        return bool(
            (self.ALDEHYDE_PATTERN is not None and mol.HasSubstructMatch(self.ALDEHYDE_PATTERN))
            or (self.KETONE_PATTERN is not None and mol.HasSubstructMatch(self.KETONE_PATTERN))
        )
