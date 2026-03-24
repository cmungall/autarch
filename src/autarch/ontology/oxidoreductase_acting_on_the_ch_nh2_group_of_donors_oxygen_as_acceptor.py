"""oxidoreductase acting on the CH-NH2 group of donors oxygen as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH2 group
acts as a hydrogen or electron donor and reduces an oxygen molecule.
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
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh2_group_of_donors import (
    OxidoreductaseActingOnTheCHNH2GroupOfDonors,
)


class OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor(
    OxidoreductaseActingOnTheCHNH2GroupOfDonors
):
    """oxidoreductase acting on the CH-NH2 group of donors oxygen as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-NH2
    group acts as a hydrogen or electron donor and reduces an oxygen molecule.
    """

    GO_ID = "GO:0016641"
    EC_NUMBER_PREFIX = "1.4.3.-"
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
    PRIMARY_AMINE_PATTERN = Chem.MolFromSmarts(
        "[N;X3,X4+;H2,H3+;!$(N=*);!$(N#*);!$(N[#6]=O);!$([N;X3,X4+]([#6])[#6])]"
    )
    AMINE_PATTERN = Chem.MolFromSmarts("[N;H1,H2,H3+;!$(N=*);!$(N#*);!$(N[#6]=O)]")
    IMINE_PATTERN = Chem.MolFromSmarts("[N;X2,X3+]=[C,N,n]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for O2-dependent CH-NH2 oxidation with peroxide formation."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH-NH2 oxidoreductase chemistry",
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

        left_reactive = self._reactive_participants(
            reaction.left_participants, self.LEFT_SPECTATOR_CHEBIS
        )
        right_reactive = self._reactive_participants(
            reaction.right_participants, self.RIGHT_SPECTATOR_CHEBIS
        )
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive participants remain after removing oxygen, peroxide, and solvent spectators",
            )

        if not any(self._is_amine_donor(participant) for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No amine-bearing donor detected among substantive substrates",
            )

        if self._has_imine_forming_oxidation(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH2 oxidoreductase: oxygen acceptor chemistry with imine-forming substrate oxidation",
            )

        if self._has_scaffold_deamination(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH2 oxidoreductase: oxygen acceptor chemistry with oxidative deamination to a carbonyl product",
            )

        if self._has_cleavage_deamination(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH2 oxidoreductase: oxygen acceptor chemistry with oxidative cleavage to a carbonyl product and a small amine coproduct",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No CH-NH2 oxidase redox signature detected",
        )

    @classmethod
    def _reactive_participants(
        cls,
        participants: list[Participant],
        spectator_chebis: set[str],
    ) -> list[Participant]:
        """Drop spectators and unstructured placeholders."""
        return [
            participant
            for participant in participants
            if participant.chebi_id not in spectator_chebis and participant.get_mol() is not None
        ]

    @classmethod
    def _has_imine_forming_oxidation(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Detect scaffold-preserving oxidation of an amine to a more oxidized nitrogen state."""
        for substrate in left_reactive:
            if not cls._has_pattern(substrate, cls.PRIMARY_AMINE_PATTERN):
                continue
            substrate_imines = cls._match_count(substrate, cls.IMINE_PATTERN)
            substrate_amines = cls._match_count(substrate, cls.AMINE_PATTERN)
            for product in right_reactive:
                if not cls._same_heavy_atom_composition(substrate, product):
                    continue
                if cls._hydrogen_count(product) >= cls._hydrogen_count(substrate):
                    continue
                if cls._match_count(product, cls.IMINE_PATTERN) > substrate_imines:
                    return True
                if cls._match_count(product, cls.AMINE_PATTERN) < substrate_amines:
                    return True
        return False

    @classmethod
    def _has_scaffold_deamination(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Detect oxidative deamination to a single carbonyl-containing product."""
        for substrate in left_reactive:
            if not cls._is_amine_donor(substrate):
                continue
            substrate_c = cls._carbon_count(substrate)
            substrate_n = cls._nitrogen_atom_count(substrate)
            substrate_o = cls._oxygen_count(substrate)
            for product in right_reactive:
                if not cls._has_pattern(product, cls.CARBONYL_PATTERN):
                    continue
                if cls._carbon_count(product) != substrate_c:
                    continue
                if cls._nitrogen_atom_count(product) != max(substrate_n - 1, 0):
                    continue
                if cls._oxygen_count(product) < substrate_o:
                    continue
                coproducts = [candidate for candidate in right_reactive if candidate is not product]
                if substrate_n == cls._nitrogen_atom_count(product):
                    continue
                if any(cls._is_small_amine_or_ammonium(candidate) for candidate in coproducts):
                    return True
        return False

    @classmethod
    def _has_cleavage_deamination(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Detect oxidative cleavage to a carbonyl fragment and a small amine fragment."""
        if len(left_reactive) != 1 or len(right_reactive) < 2:
            return False

        substrate = left_reactive[0]
        if not cls._is_amine_donor(substrate):
            return False

        substrate_c = cls._carbon_count(substrate)
        for carbonyl_product in right_reactive:
            if not cls._has_pattern(carbonyl_product, cls.CARBONYL_PATTERN):
                continue
            for nitrogen_product in right_reactive:
                if nitrogen_product is carbonyl_product:
                    continue
                if not cls._is_small_amine_or_ammonium(nitrogen_product):
                    continue
                if (
                    cls._carbon_count(carbonyl_product) + cls._carbon_count(nitrogen_product)
                    == substrate_c
                ):
                    return True
        return False

    @classmethod
    def _is_amine_donor(cls, participant: Participant) -> bool:
        return cls._has_pattern(participant, cls.AMINE_PATTERN)

    @classmethod
    def _is_small_amine_or_ammonium(cls, participant: Participant) -> bool:
        if participant.chebi_id in {CHEBI_NH3, CHEBI_NH4}:
            return True
        return (
            cls._nitrogen_atom_count(participant) > 0
            and cls._carbon_count(participant) <= 2
            and not cls._has_pattern(participant, cls.CARBONYL_PATTERN)
        )

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

    @classmethod
    def _same_heavy_atom_composition(cls, left: Participant, right: Participant) -> bool:
        """Compare heavy-atom composition."""
        return cls._heavy_atom_counts(left) == cls._heavy_atom_counts(right)

    @staticmethod
    def _hydrogen_count(participant: Participant) -> int:
        """Count represented hydrogens."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        """Count carbons."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_atom_count(participant: Participant) -> int:
        """Count nitrogens."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        """Count oxygens."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @staticmethod
    def _match_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count SMARTS matches in a participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))

    @classmethod
    def _has_pattern(cls, participant: Participant, pattern: Chem.Mol | None) -> bool:
        """Check whether a participant matches a SMARTS pattern."""
        return cls._match_count(participant, pattern) > 0
