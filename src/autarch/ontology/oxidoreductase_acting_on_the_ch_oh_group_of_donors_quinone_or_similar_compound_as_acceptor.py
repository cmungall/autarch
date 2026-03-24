"""oxidoreductase acting on the CH-OH group of donors, quinone or similar compound as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-OH group acts as a hydrogen or electron donor and reduces a quinone or a similar acceptor molecule.
"""

from __future__ import annotations

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
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE


class OxidoreductaseActingOnTheCHOHGroupOfDonorsQuinoneOrSimilarCompoundAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase acting on the CH-OH group of donors, quinone or similar compound as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-OH group acts as a hydrogen or electron donor and reduces a quinone or a similar acceptor molecule.
    """

    GO_ID = "GO:0016901"
    EC_NUMBER_PREFIX = "1.1.5.-"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    QUINONE_ACCEPTOR_IDS = {"CHEBI:132124"}
    QUINOL_PRODUCT_IDS = {"CHEBI:24646"}
    SPECTATOR_IDS = {CHEBI_H_PLUS, CHEBI_H2O, CHEBI_H2O2}
    NICOTINAMIDE_IDS = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1,-][C;!$(C=O)]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CH-OH/quinone redox chemistry",
            )

        all_ids = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if CHEBI_O2 in all_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen-dependent chemistry indicates a different oxidoreductase branch",
            )
        if self.NICOTINAMIDE_IDS & all_ids:
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P)-linked chemistry belongs to a different CH-OH oxidoreductase branch",
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
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        if not any(participant.chebi_id in self.QUINONE_ACCEPTOR_IDS for participant in left):
            return ClassificationResult(
                is_member=False,
                explanation="Requires a quinone or similar acceptor on the donor side",
            )
        if not any(participant.chebi_id in self.QUINOL_PRODUCT_IDS for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="Requires the reduced quinol-like acceptor product",
            )

        left_core = self._reactive_participants(left, self.QUINONE_ACCEPTOR_IDS)
        right_core = self._reactive_participants(right, self.QUINOL_PRODUCT_IDS)
        donor = self._collapse_uniform_participants(left_core)
        product = self._collapse_uniform_participants(right_core)
        if donor is None or product is None:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one substantive CH-OH donor and one oxidized product",
            )
        if donor.get_mol() is None or product.get_mol() is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for CH-OH donor/product comparison",
            )
        if not self._same_heavy_atom_composition(donor, product):
            return ClassificationResult(
                is_member=False,
                explanation="Substantive donor/product pair does not preserve the same heavy-atom scaffold",
            )
        if self._pattern_count(donor, self.HYDROXYL_PATTERN) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Substantive donor does not contain a supported CH-OH group",
            )
        if self._pattern_count(product, self.CARBONYL_PATTERN) <= self._pattern_count(
            donor, self.CARBONYL_PATTERN
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Substantive product does not increase carbonyl content",
            )
        if self._hydrogen_count(product) >= self._hydrogen_count(donor):
            return ClassificationResult(
                is_member=False,
                explanation="Substantive product is not more oxidized than the donor",
            )

        return ClassificationResult(
            is_member=True,
            explanation="CH-OH oxidoreductase (quinone branch): oxidation of a hydroxyl-bearing donor coupled to quinone reduction",
        )

    def _reactive_participants(
        self,
        participants: list[Participant],
        excluded_ids: set[str],
    ) -> list[Participant]:
        return [
            participant
            for participant in participants
            if participant.chebi_id not in excluded_ids
            if participant.chebi_id not in self.SPECTATOR_IDS
            if participant.chebi_id is not None or participant.smiles is not None
        ]

    @staticmethod
    def _collapse_uniform_participants(participants: list[Participant]) -> Participant | None:
        if not participants:
            return None
        representative = participants[0]
        if all(representative.is_same_molecule(candidate) for candidate in participants[1:]):
            return representative
        return None

    @staticmethod
    def _heavy_atom_counts(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() > 1
        )

    @classmethod
    def _same_heavy_atom_composition(cls, left: Participant, right: Participant) -> bool:
        return cls._heavy_atom_counts(left) == cls._heavy_atom_counts(right)

    @staticmethod
    def _hydrogen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(atom.GetTotalNumHs(includeNeighbors=True) for atom in mol.GetAtoms())

    @staticmethod
    def _pattern_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))
