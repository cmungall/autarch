"""glycine oxidase activity.

Catalysis of oxidation of glycine and related small amino acids with dioxygen as
acceptor, yielding an oxo acid, hydrogen peroxide, and a small amine product.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import CHEBI_H2O, CHEBI_H2O2, CHEBI_H_PLUS, CHEBI_NH4, CHEBI_O2
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh2_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor,
)

CHEBI_GLYOXYLATE = "CHEBI:36655"
CHEBI_PYRUVATE = "CHEBI:15361"


class GlycineOxidase(OxidoreductaseActingOnTheCHNH2GroupOfDonorsOxygenAsAcceptor):
    """glycine oxidase activity.

    Catalysis of oxidation of glycine and related small amino acids with
    dioxygen as acceptor, yielding an oxo acid, hydrogen peroxide, and a small
    amine product.
    """

    GO_ID = "GO:0043799"
    EC_NUMBER_PREFIX = "1.4.3.19"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for oxidation of a small amino acid substrate to glyoxylate or pyruvate."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        if not any(participant.chebi_id == CHEBI_O2 for participant in reaction.left_participants):
            return ClassificationResult(is_member=False, explanation="No dioxygen reactant found")
        if not any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return ClassificationResult(is_member=False, explanation="No water reactant found")
        if not any(participant.chebi_id == CHEBI_H2O2 for participant in reaction.right_participants):
            return ClassificationResult(is_member=False, explanation="No hydrogen peroxide product found")

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in {CHEBI_O2, CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in {CHEBI_H2O2, CHEBI_H_PLUS}
        ]

        if not any(self._is_small_amino_acid_substrate(participant) for participant in left_core):
            return ClassificationResult(
                is_member=False,
                explanation="No small amino-acid substrate detected",
            )

        has_oxo_acid = any(self._is_allowed_oxo_acid(participant) for participant in right_core)
        if not has_oxo_acid:
            return ClassificationResult(
                is_member=False,
                explanation="No glyoxylate- or pyruvate-like product detected",
            )

        has_small_amine = any(self._is_small_amine_product(participant) for participant in right_core)
        if not has_small_amine:
            return ClassificationResult(
                is_member=False,
                explanation="No ammonium or small amine coproduct detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Glycine oxidase: oxidation of a small amino acid to glyoxylate/pyruvate with peroxide formation",
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @classmethod
    def _is_small_amino_acid_substrate(cls, participant: Participant) -> bool:
        return bool(
            participant.has_moiety(Moiety.CARBOXYL)
            and cls._nitrogen_count(participant) >= 1
            and 2 <= cls._carbon_count(participant) <= 4
            and not participant.is_thioester()
            and not participant.is_phosphorylated()
        )

    @classmethod
    def _is_allowed_oxo_acid(cls, participant: Participant) -> bool:
        if participant.chebi_id in {CHEBI_GLYOXYLATE, CHEBI_PYRUVATE}:
            return True
        return bool(
            participant.has_moiety(Moiety.CARBOXYL)
            and cls._carbon_count(participant) <= 3
            and cls._nitrogen_count(participant) == 0
        )

    @classmethod
    def _is_small_amine_product(cls, participant: Participant) -> bool:
        if participant.chebi_id == CHEBI_NH4:
            return True
        return bool(cls._nitrogen_count(participant) >= 1 and cls._carbon_count(participant) <= 2)
