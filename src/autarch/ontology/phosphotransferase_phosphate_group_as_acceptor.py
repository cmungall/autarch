"""phosphotransferase phosphate group as acceptor.

Catalysis of the transfer of a phosphorus-containing group from one compound
(donor) to a phosphate group (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CTP,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_UDP,
    CHEBI_UTP,
)
from autarch.ontology.transferase import Transferase


class PhosphotransferasePhosphateGroupAsAcceptor(Transferase):
    """phosphotransferase phosphate group as acceptor.

    Catalysis of the transfer of a phosphorus-containing group from one
    compound (donor) to a phosphate group (acceptor).
    """

    GO_ID = "GO:0016776"
    EC_NUMBER_PREFIX = "2.7.4.-"
    DONOR_PRODUCT_PAIRS = {
        CHEBI_ATP: CHEBI_ADP,
        CHEBI_GTP: CHEBI_GDP,
        CHEBI_UTP: CHEBI_UDP,
        CHEBI_CTP: CHEBI_CDP,
    }
    EXCLUDED_CHEBIS = {
        CHEBI_H2O,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_O2,
    }
    PHOSPHORIC_ANHYDRIDE_PATTERN = Chem.MolFromSmarts("[P][O][P]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for phosphoryl transfer onto a pre-existing phosphate group."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not phosphate-acceptor phosphotransferase chemistry",
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
        """Evaluate one phosphotransfer direction."""
        chebis = {
            participant.chebi_id
            for participant in left + right
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Hydrolysis, oxygenation, or redox cofactors indicate different chemistry",
            )

        donor_pair = self._find_donor_pair(left, right)
        if donor_pair is None:
            return ClassificationResult(
                is_member=False,
                explanation="No nucleoside triphosphate to diphosphate donor pair found",
            )

        donor_id, donor_product_id = donor_pair
        left_reactive = self._drop_one_chebi(left, donor_id)
        right_reactive = self._drop_one_chebi(right, donor_product_id)
        left_reactive = [
            participant for participant in left_reactive if participant.chebi_id != CHEBI_H_PLUS
        ]
        right_reactive = [
            participant for participant in right_reactive if participant.chebi_id != CHEBI_H_PLUS
        ]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive phosphate acceptor remains after removing donor nucleotides",
            )

        if any(participant.get_mol() is None for participant in left_reactive + right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Missing structural information for phosphate acceptor substrates or products",
            )

        left_phosphate_acceptors = [
            participant for participant in left_reactive if self._phosphorus_count(participant) >= 1
        ]
        if not left_phosphate_acceptors:
            return ClassificationResult(
                is_member=False,
                explanation="Acceptor substrate lacks a pre-existing phosphate group",
            )

        left_total_p = sum(self._phosphorus_count(participant) for participant in left_reactive)
        right_total_p = sum(self._phosphorus_count(participant) for participant in right_reactive)
        if right_total_p != left_total_p + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Non-donor participants do not gain exactly one phosphorus atom",
            )

        left_anhydrides = sum(self._phosphoric_anhydride_count(participant) for participant in left_reactive)
        right_anhydrides = sum(self._phosphoric_anhydride_count(participant) for participant in right_reactive)
        if right_anhydrides <= left_anhydrides:
            return ClassificationResult(
                is_member=False,
                explanation="Acceptor products do not form a new phosphoric anhydride bond",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Phosphotransferase with phosphate acceptor: a nucleoside triphosphate donates one phosphoryl group to an already phosphorylated acceptor",
        )

    @classmethod
    def _find_donor_pair(
        cls,
        left: list[Participant],
        right: list[Participant],
    ) -> tuple[str, str] | None:
        left_ids = {
            participant.chebi_id
            for participant in left
            if participant.chebi_id
        }
        right_ids = {
            participant.chebi_id
            for participant in right
            if participant.chebi_id
        }
        for donor_id, donor_product_id in cls.DONOR_PRODUCT_PAIRS.items():
            if donor_id in left_ids and donor_product_id in right_ids:
                return donor_id, donor_product_id
        return None

    @staticmethod
    def _drop_one_chebi(participants: list[Participant], chebi_id: str) -> list[Participant]:
        remaining: list[Participant] = []
        removed = False
        for participant in participants:
            if not removed and participant.chebi_id == chebi_id:
                if participant.count > 1:
                    remaining.append(participant.model_copy(update={"count": participant.count - 1}))
                removed = True
                continue
            remaining.append(participant)
        return remaining

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _phosphoric_anhydride_count(cls, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None or cls.PHOSPHORIC_ANHYDRIDE_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(cls.PHOSPHORIC_ANHYDRIDE_PATTERN))
