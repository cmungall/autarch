"""oxidoreductase acting on the CH-NH group of donors NAD or NADP as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH group acts
as a hydrogen or electron donor and reduces NAD or NADP.
"""

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on the CH-NH group of donors NAD or NADP as
    acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH group
    acts as a hydrogen or electron donor and reduces NAD or NADP.
    """

    GO_ID = "GO:0016646"
    EC_NUMBER_PREFIX = "1.5.1.-"
    OXIDIZED_NICOTINAMIDES = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    REDUCED_NICOTINAMIDES = {CHEBI_NADH, CHEBI_NADPH}
    LEFT_SPECTATOR_CHEBIS = OXIDIZED_NICOTINAMIDES | {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = REDUCED_NICOTINAMIDES | {CHEBI_H2O, CHEBI_H_PLUS}
    SECONDARY_AMINE_PATTERN = Chem.MolFromSmarts(
        "[N;X3,X4+;H0,H1,H2;!$(N[#6]=O)]([#6])[#6]"
    )
    AMINE_PATTERN = Chem.MolFromSmarts("[N;H1,H2,H3+;!$(N=*);!$(N#*)]")
    IMINE_PATTERN = Chem.MolFromSmarts("[N;X2,X3+]=[C,N,n]")
    AROMATIC_N_PATTERN = Chem.MolFromSmarts("[nH0]")
    AROMATIC_NH_PATTERN = Chem.MolFromSmarts("[nH1]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for NAD(P)-linked CH-NH oxidoreduction."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH-NH oxidoreductase",
            )

        if self._has_external_oxidant(reaction):
            return ClassificationResult(
                is_member=False,
                explanation="External O2/H2O2 oxidant indicates a different oxidoreductase branch",
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
                explanation="NAD(P)-linked CH-NH oxidoreduction detected in the reverse direction",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No NAD(P)-linked CH-NH redox signature detected",
        )

    def _check_direction(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        """Check one oxidation direction: reduced substrate plus NAD(P)+."""
        if not self._has_nicotinamide_pair(left, right):
            return ClassificationResult(
                is_member=False,
                explanation="Requires NAD(+) or NADP(+) reduced across the reaction",
            )

        left_reactive = self._reactive_participants(left, self.LEFT_SPECTATOR_CHEBIS)
        right_reactive = self._reactive_participants(right, self.RIGHT_SPECTATOR_CHEBIS)
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive non-cofactor participants remain after removing spectators",
            )

        if any(participant.get_mol() is None for participant in left_reactive + right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for substantive participants",
            )

        if any(participant.chebi_id in {CHEBI_NH3, CHEBI_NH4} for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Free ammonia production indicates CH-NH2 dehydrogenase chemistry",
            )

        if self._has_scaffold_preserving_n_redox(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)-linked CH-NH oxidoreduction with a scaffold-preserving nitrogen redox transformation",
            )

        if self._has_hydrolytic_chnh_cleavage(left, left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)-linked CH-NH oxidoreduction with water-assisted cleavage of a secondary amine substrate",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No scaffold-preserving or cleavage CH-NH redox mode matched",
        )

    def _has_scaffold_preserving_n_redox(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Detect nitrogen-centered redox without substrate fragmentation."""
        for substrate in left_reactive:
            for product in right_reactive:
                if not self._same_heavy_atom_composition(substrate, product):
                    continue
                if self._nitrogen_atom_count(substrate) == 0:
                    continue
                hydrogen_loss = self._hydrogen_count(substrate) - self._hydrogen_count(product)
                if hydrogen_loss <= 0 or hydrogen_loss > 4:
                    continue
                if self._nitrogen_redox_progression(substrate, product):
                    return True
        return False

    def _has_hydrolytic_chnh_cleavage(
        self,
        left: list[Participant],
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Detect opine/saccharopine-like cleavage of a secondary amine substrate."""
        if not any(participant.chebi_id == CHEBI_H2O for participant in left):
            return False
        if len(left_reactive) != 1 or len(right_reactive) < 2:
            return False

        substrate = left_reactive[0]
        if not self._has_pattern(substrate, self.SECONDARY_AMINE_PATTERN):
            return False

        has_n_product = any(self._nitrogen_atom_count(participant) > 0 for participant in right_reactive)
        has_carbonyl_product = any(
            self._has_pattern(participant, self.CARBONYL_PATTERN) for participant in right_reactive
        )
        return has_n_product and has_carbonyl_product

    def _has_nicotinamide_pair(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        """Require a clean NAD(P)+ -> NAD(P)H conversion."""
        left_oxidized = self._count_chebis(left, self.OXIDIZED_NICOTINAMIDES)
        right_reduced = self._count_chebis(right, self.REDUCED_NICOTINAMIDES)
        left_reduced = self._count_chebis(left, self.REDUCED_NICOTINAMIDES)
        right_oxidized = self._count_chebis(right, self.OXIDIZED_NICOTINAMIDES)
        return left_oxidized in {1, 2} and left_oxidized == right_reduced and left_reduced == 0 and right_oxidized == 0

    @staticmethod
    def _has_external_oxidant(reaction: Reaction) -> bool:
        """Reject oxygen/peroxide-driven oxidase or oxygenase chemistry."""
        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        return CHEBI_O2 in chebis or CHEBI_H2O2 in chebis

    @staticmethod
    def _reactive_participants(
        participants: list[Participant],
        spectator_chebis: set[str],
    ) -> list[Participant]:
        """Drop nicotinamides, water, and hydrons from the reactive core."""
        return [
            participant
            for participant in participants
            if participant.chebi_id not in spectator_chebis
        ]

    @staticmethod
    def _count_chebis(participants: list[Participant], chebis: set[str]) -> int:
        """Count a set of ChEBI IDs with stoichiometry."""
        return sum(
            max(participant.count, 1)
            for participant in participants
            if participant.chebi_id in chebis
        )

    @classmethod
    def _nitrogen_redox_progression(
        cls,
        substrate: Participant,
        product: Participant,
    ) -> bool:
        """Check that nitrogen moves from a more reduced to a more oxidized state."""
        substrate_profile = cls._nitrogen_profile(substrate)
        product_profile = cls._nitrogen_profile(product)
        substrate_reduced = substrate_profile["amine"] + substrate_profile["aromatic_nh"]
        product_reduced = product_profile["amine"] + product_profile["aromatic_nh"]
        substrate_oxidized = substrate_profile["imine"] + substrate_profile["aromatic_n"]
        product_oxidized = product_profile["imine"] + product_profile["aromatic_n"]
        if substrate_reduced > product_reduced and product_oxidized > substrate_oxidized:
            return True
        return (
            substrate_reduced - product_reduced >= 2
            and product_oxidized >= substrate_oxidized
        )

    @classmethod
    def _nitrogen_profile(cls, participant: Participant) -> Counter[str]:
        """Summarize nitrogen environments relevant to CH-NH redox."""
        return Counter(
            amine=cls._match_count(participant, cls.AMINE_PATTERN),
            imine=cls._match_count(participant, cls.IMINE_PATTERN),
            aromatic_n=cls._match_count(participant, cls.AROMATIC_N_PATTERN),
            aromatic_nh=cls._match_count(participant, cls.AROMATIC_NH_PATTERN),
        )

    @classmethod
    def _same_heavy_atom_composition(cls, left: Participant, right: Participant) -> bool:
        """Compare atom counts excluding hydrogen."""
        return cls._heavy_atom_counts(left) == cls._heavy_atom_counts(right)

    @staticmethod
    def _heavy_atom_counts(participant: Participant) -> Counter[int]:
        """Count heavy atoms by atomic number."""
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() != 1
        )

    @staticmethod
    def _hydrogen_count(participant: Participant) -> int:
        """Count total hydrogens represented on the molecule."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

    @staticmethod
    def _nitrogen_atom_count(participant: Participant) -> int:
        """Count nitrogen atoms in a participant."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _match_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count SMARTS matches for a single participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))

    @classmethod
    def _has_pattern(cls, participant: Participant, pattern: Chem.Mol | None) -> bool:
        """Check if a participant has a SMARTS pattern."""
        return cls._match_count(participant, pattern) > 0
