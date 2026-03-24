"""oxidoreductase acting on NAD(P)H, quinone or similar compound as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which NADH or NADPH
acts as a hydrogen or electron donor and reduces a quinone or a similar
acceptor molecule.
"""

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on NAD(P)H, quinone or similar compound as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which NADH or
    NADPH acts as a hydrogen or electron donor and reduces a quinone or a
    similar acceptor molecule.
    """

    GO_ID = "GO:0016655"
    EC_NUMBER_PREFIX = "1.6.5.-"
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    SPECIAL_ACCEPTOR_IDS = {"CHEBI:59513"}  # monodehydro-L-ascorbate
    SPECIAL_REDUCED_PRODUCT_IDS = {"CHEBI:38290"}  # L-ascorbate
    HYDROXYL_PATTERN = Chem.MolFromSmarts("[O;H1,-][C;!$(C=O)]")
    CARBONYL_PATTERN = Chem.MolFromSmarts("[CX3]=[OX1]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for quinone-like reduction by NAD(P)H."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not quinone reduction chemistry",
            )

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        return forward

    def _check_direction(
        self,
        left: list[Participant],
        right: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        if any(participant.chebi_id == CHEBI_O2 for participant in left + right):
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen-dependent chemistry indicates a different oxidoreductase branch",
            )

        donor_pair = self._find_cofactor_pair(left, right)
        if donor_pair is None:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD(P)H to NAD(P)+ donor/product pair found",
            )

        donor_id, oxidized_id = donor_pair
        left_reactive = self._drop_one_chebi(left, donor_id)
        right_reactive = self._drop_one_chebi(right, oxidized_id)
        left_reactive = [
            participant
            for participant in left_reactive
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_reactive
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        substrate = self._collapse_uniform_participants(left_reactive)
        product = self._collapse_uniform_participants(right_reactive)
        if substrate is None or product is None:
            return ClassificationResult(
                is_member=False,
                explanation="Quinone reduction requires one substantive acceptor and one substantive reduced product",
            )
        if substrate.get_mol() is None or product.get_mol() is None:
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for quinone-like substrates or products",
            )

        if not self._same_heavy_atom_composition(substrate, product):
            return ClassificationResult(
                is_member=False,
                explanation="Reactive substrate and product do not preserve the same heavy-atom scaffold",
            )

        if self._hydrogen_count(product) <= self._hydrogen_count(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Reactive product is not more reduced than the acceptor substrate",
            )

        if self._is_special_similar_acceptor_pair(substrate, product):
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)H-dependent reduction of a similar acceptor scaffold",
            )

        if not self._is_quinone_like_acceptor(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Substantive substrate is not quinone-like",
            )

        substrate_oxygen = self._oxygen_count(substrate)
        product_oxygen = self._oxygen_count(product)
        if substrate_oxygen != product_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="Quinone reduction should preserve the oxygen inventory of the acceptor scaffold",
            )

        substrate_carbonyls = self._pattern_count(substrate, self.CARBONYL_PATTERN)
        product_carbonyls = self._pattern_count(product, self.CARBONYL_PATTERN)
        if product_carbonyls >= substrate_carbonyls:
            return ClassificationResult(
                is_member=False,
                explanation="Reduced product does not decrease quinone-like carbonyl content",
            )

        return ClassificationResult(
            is_member=True,
            explanation="NAD(P)H-dependent reduction of a quinone-like acceptor scaffold",
        )

    @staticmethod
    def _find_cofactor_pair(
        left: list[Participant],
        right: list[Participant],
    ) -> tuple[str, str] | None:
        """Find the nicotinamide donor/product pair for the orientation."""
        left_ids = {participant.chebi_id for participant in left if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right if participant.chebi_id}
        pairs = (
            (CHEBI_NADH, CHEBI_NAD_PLUS),
            (CHEBI_NADPH, CHEBI_NADP_PLUS),
        )
        for donor_id, oxidized_id in pairs:
            if donor_id in left_ids and oxidized_id in right_ids:
                return donor_id, oxidized_id
        return None

    @staticmethod
    def _drop_one_chebi(participants: list[Participant], chebi_id: str) -> list[Participant]:
        """Remove one instance of a participant, respecting stoichiometric count."""
        remaining: list[Participant] = []
        removed = False
        for participant in participants:
            if not removed and participant.chebi_id == chebi_id:
                if participant.count > 1:
                    remaining.append(
                        participant.model_copy(update={"count": participant.count - 1})
                    )
                removed = True
                continue
            remaining.append(participant)
        return remaining

    @staticmethod
    def _collapse_uniform_participants(participants: list[Participant]) -> Participant | None:
        """Collapse repeated copies of the same participant into one representative."""
        if not participants:
            return None
        representative = participants[0]
        if all(representative.is_same_molecule(candidate) for candidate in participants[1:]):
            return representative
        return None

    def _is_special_similar_acceptor_pair(
        self,
        substrate: Participant,
        product: Participant,
    ) -> bool:
        """Handle similar non-quinone acceptors represented cleanly in RHEA."""
        return (
            substrate.chebi_id in self.SPECIAL_ACCEPTOR_IDS
            and product.chebi_id in self.SPECIAL_REDUCED_PRODUCT_IDS
        )

    def _is_quinone_like_acceptor(self, participant: Participant) -> bool:
        """Detect quinone-like oxidized acceptors by ring and carbonyl content."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return (
            mol.GetRingInfo().NumRings() >= 1
            and self._oxygen_count(participant) >= 2
            and self._pattern_count(participant, self.CARBONYL_PATTERN) >= 2
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

    def _same_heavy_atom_composition(self, left: Participant, right: Participant) -> bool:
        """Compare heavy-atom composition between two participants."""
        return self._heavy_atom_counts(left) == self._heavy_atom_counts(right)

    @staticmethod
    def _hydrogen_count(participant: Participant) -> int:
        """Count explicit and implicit hydrogens."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

    @staticmethod
    def _oxygen_count(participant: Participant) -> int:
        """Count oxygen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    @staticmethod
    def _pattern_count(participant: Participant, pattern: Chem.Mol | None) -> int:
        """Count substructure matches for a participant."""
        mol = participant.get_mol()
        if mol is None or pattern is None:
            return 0
        return len(mol.GetSubstructMatches(pattern))
