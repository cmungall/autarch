"""oxidoreductase acting on the CH-NH group of donors oxygen as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH group
acts as a hydrogen or electron donor and reduces oxygen.
"""

from collections import Counter

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
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh_group_of_donors import (
    OxidoreductaseActingOnTheCHNHGroupOfDonors,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor as CHNHNADPBase,
)


class OxidoreductaseActingOnTheCHNHGroupOfDonorsOxygenAsAcceptor(
    OxidoreductaseActingOnTheCHNHGroupOfDonors
):
    """oxidoreductase acting on the CH-NH group of donors oxygen as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH
    group acts as a hydrogen or electron donor and reduces oxygen.
    """

    GO_ID = "GO:0016647"
    EC_NUMBER_PREFIX = "1.5.3.-"
    SECONDARY_AMINE_PATTERN = CHNHNADPBase.SECONDARY_AMINE_PATTERN
    TERTIARY_AMINE_PATTERN = Chem.MolFromSmarts(
        "[N;X3,X4+;H0,H1;!$(N[#6]=O)]([#6])([#6])[#6]"
    )
    CARBONYL_PATTERN = CHNHNADPBase.CARBONYL_PATTERN
    ALPHA_AMINO_ACID_PATTERN = Chem.MolFromSmarts("[N;X3,X4+][CH1,CH2][C](=O)[O;H1,X1-]")
    OXO_ACID_PATTERN = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)[O;H1,X1-]")
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O2, CHEBI_H2O, CHEBI_H_PLUS}
    EXCLUDED_SOLUBLE_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for oxygen-dependent CH-NH oxidation with peroxide formation."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH-NH oxidase chemistry",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen reactant found",
            )

        if not any(participant.chebi_id == CHEBI_H2O2 for participant in reaction.right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No hydrogen peroxide product found",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_SOLUBLE_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Soluble nicotinamide or flavin cofactors indicate a different oxidoreductase branch",
            )

        left_reactive = CHNHNADPBase._reactive_participants(
            reaction.left_participants, self.LEFT_SPECTATOR_CHEBIS
        )
        right_reactive = CHNHNADPBase._reactive_participants(
            reaction.right_participants, self.RIGHT_SPECTATOR_CHEBIS
        )
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive CH-NH donor or oxidized products remain after removing spectators",
            )

        if any(participant.get_mol() is None for participant in left_reactive + right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for CH-NH oxidase chemistry",
            )

        if not any(self._is_chnh_donor(participant) for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No CH-NH donor detected among substantive substrates",
            )

        if self._is_amino_acid_oxidative_deamination(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Amino-acid oxidative deamination belongs to the CH-NH2 oxidase branch",
            )

        if self._has_scaffold_preserving_n_redox(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH oxidase: oxygen acceptor chemistry with hydrogen peroxide formation from a nitrogen-bearing donor",
            )

        if self._has_hydrolytic_chnh_cleavage(reaction.left_participants, left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH oxidase: oxygen acceptor chemistry with oxidative cleavage of a CH-NH substrate and hydrogen peroxide formation",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No CH-NH oxidase redox signature detected",
        )

    @classmethod
    def _is_chnh_donor(cls, participant: Participant) -> bool:
        return CHNHNADPBase._has_pattern(participant, cls.SECONDARY_AMINE_PATTERN) or CHNHNADPBase._has_pattern(
            participant, cls.TERTIARY_AMINE_PATTERN
        )

    def _has_scaffold_preserving_n_redox(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
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
        if not any(participant.chebi_id == CHEBI_H2O for participant in left):
            return False
        if len(left_reactive) != 1 or not right_reactive:
            return False

        substrate = left_reactive[0]
        if not self._is_chnh_donor(substrate):
            return False

        has_n_product = any(self._nitrogen_atom_count(participant) > 0 for participant in right_reactive)
        has_carbonyl_product = any(
            self._has_pattern(participant, self.CARBONYL_PATTERN) for participant in right_reactive
        )
        return has_n_product and has_carbonyl_product

    @classmethod
    def _nitrogen_redox_progression(cls, substrate: Participant, product: Participant) -> bool:
        substrate_profile = cls._nitrogen_profile(substrate)
        product_profile = cls._nitrogen_profile(product)
        substrate_reduced = (
            substrate_profile["amine"]
            + substrate_profile["tertiary_amine"]
            + substrate_profile["aromatic_nh"]
        )
        product_reduced = (
            product_profile["amine"]
            + product_profile["tertiary_amine"]
            + product_profile["aromatic_nh"]
        )
        substrate_oxidized = substrate_profile["imine"] + substrate_profile["aromatic_n"]
        product_oxidized = product_profile["imine"] + product_profile["aromatic_n"]
        if substrate_reduced > product_reduced and product_oxidized > substrate_oxidized:
            return True
        return (
            substrate_reduced - product_reduced >= 1
            and product_oxidized > substrate_oxidized
        )

    @classmethod
    def _nitrogen_profile(cls, participant: Participant):
        return Counter(
            amine=cls._match_count(participant, CHNHNADPBase.AMINE_PATTERN),
            tertiary_amine=cls._match_count(participant, cls.TERTIARY_AMINE_PATTERN),
            imine=cls._match_count(participant, CHNHNADPBase.IMINE_PATTERN),
            aromatic_n=cls._match_count(participant, CHNHNADPBase.AROMATIC_N_PATTERN),
            aromatic_nh=cls._match_count(participant, CHNHNADPBase.AROMATIC_NH_PATTERN),
        )

    def _is_amino_acid_oxidative_deamination(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        return (
            len(left_reactive) == 1
            and self._has_pattern(left_reactive[0], self.ALPHA_AMINO_ACID_PATTERN)
            and any(self._has_pattern(participant, self.OXO_ACID_PATTERN) for participant in right_reactive)
            and any(self._is_small_amine(participant) for participant in right_reactive)
        )

    @staticmethod
    def _same_heavy_atom_composition(left: Participant, right: Participant) -> bool:
        return CHNHNADPBase._same_heavy_atom_composition(left, right)

    @staticmethod
    def _heavy_atom_counts(participant: Participant):
        return CHNHNADPBase._heavy_atom_counts(participant)

    @staticmethod
    def _hydrogen_count(participant: Participant) -> int:
        return CHNHNADPBase._hydrogen_count(participant)

    @staticmethod
    def _nitrogen_atom_count(participant: Participant) -> int:
        return CHNHNADPBase._nitrogen_atom_count(participant)

    @staticmethod
    def _match_count(participant: Participant, pattern) -> int:
        return CHNHNADPBase._match_count(participant, pattern)

    @classmethod
    def _has_pattern(cls, participant: Participant, pattern) -> bool:
        return CHNHNADPBase._has_pattern(participant, pattern)

    @staticmethod
    def _is_small_amine(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        return nitrogen_count == 1 and oxygen_count == 0 and carbon_count <= 2
