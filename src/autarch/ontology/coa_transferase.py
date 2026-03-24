"""CoA-transferase.

Catalysis of the transfer of a coenzyme A (CoA) group from one compound
(donor) to another (acceptor).
"""

from collections import Counter

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_GTP,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.moiety import Moiety
from autarch.ontology.transferase import Transferase


class CoATransferase(Transferase):
    """CoA-transferase.

    Catalysis of the transfer of a coenzyme A (CoA) group from one compound
    (donor) to another (acceptor).
    """

    GO_ID = "GO:0008410"
    EC_NUMBER_PREFIX = "2.8.3.-"
    EXCLUDED_CHEBIS = {
        CHEBI_ATP,
        CHEBI_GTP,
        CHEBI_COA,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_O2,
        CHEBI_SAH,
        CHEBI_SAM,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for CoA transfer between a donor thioester and an acceptor acid."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not CoA transferase chemistry",
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
        """Evaluate one CoA transfer direction."""
        chebis = {
            participant.chebi_id
            for participant in left + right
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Water, free CoA, nucleotide, or redox cofactors indicate different chemistry",
            )

        if any(participant.get_mol() is None for participant in left + right):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for CoA-transferase substrates or products",
            )

        left_thioesters = [participant for participant in left if self._is_coa_like_thioester(participant)]
        right_thioesters = [participant for participant in right if self._is_coa_like_thioester(participant)]
        if len(left_thioesters) != 1 or len(right_thioesters) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="CoA transfer requires one CoA-thioester donor and one CoA-thioester product",
            )

        left_acids = [participant for participant in left if self._is_non_thioester_carboxylate(participant)]
        right_acids = [participant for participant in right if self._is_non_thioester_carboxylate(participant)]
        if not left_acids or not right_acids:
            return ClassificationResult(
                is_member=False,
                explanation="CoA transfer requires non-thioester carboxylate acceptor and donor-acid products",
            )

        if self._heavy_atom_counts(left_thioesters[0]) == self._heavy_atom_counts(right_thioesters[0]):
            return ClassificationResult(
                is_member=False,
                explanation="Thioester scaffold does not change across the reaction",
            )

        return ClassificationResult(
            is_member=True,
            explanation="CoA-transferase: CoA-thioester donor exchanges its acyl group with a carboxylate acceptor",
        )

    @classmethod
    def _is_coa_like_thioester(cls, participant: Participant) -> bool:
        return (
            participant.is_thioester()
            and cls._atom_count(participant, 15) >= 3
            and cls._atom_count(participant, 7) >= 3
        )

    @staticmethod
    def _is_non_thioester_carboxylate(participant: Participant) -> bool:
        return participant.has_moiety(Moiety.CARBOXYL) and not participant.is_thioester()

    @staticmethod
    def _atom_count(participant: Participant, atomic_num: int) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)

    @staticmethod
    def _heavy_atom_counts(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() != 1
        )
