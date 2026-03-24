"""phosphotransferase activity, nitrogenous group as acceptor.

Catalysis of the transfer of a phosphorus-containing group from one compound
(donor) to a nitrogenous group (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
)
from autarch.ontology.transferase import Transferase

CHEBI_PEP = "CHEBI:58702"
CHEBI_PYRUVATE = "CHEBI:15361"


class PhosphotransferaseNitrogenousGroupAsAcceptor(Transferase):
    """phosphotransferase activity, nitrogenous group as acceptor.

    Catalysis of the transfer of a phosphorus-containing group from one
    compound (donor) to a nitrogenous group (acceptor).
    """

    GO_ID = "GO:0016775"
    EC_NUMBER_PREFIX = "2.7.3.-"
    DONOR_PRODUCT_PAIRS = {
        CHEBI_ATP: CHEBI_ADP,
        CHEBI_GTP: CHEBI_GDP,
        CHEBI_PEP: CHEBI_PYRUVATE,
    }
    LEFT_SPECTATOR_CHEBIS = set(DONOR_PRODUCT_PAIRS)
    RIGHT_SPECTATOR_CHEBIS = set(DONOR_PRODUCT_PAIRS.values()) | {CHEBI_H2O, CHEBI_H_PLUS}
    PHOSPHORAMIDATE_PATTERN = Chem.MolFromSmarts("[n,N][P](=[O])([O])[O]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for phosphoryl transfer onto a nitrogenous group."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not phosphotransfer to a nitrogenous group",
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
        donor_pair = self._find_donor_pair(left, right)
        if donor_pair is None:
            return ClassificationResult(
                is_member=False,
                explanation="No phosphorus-donor to spent-donor pair found",
            )

        donor_id, donor_product_id = donor_pair
        left_reactive = self._drop_one_chebi(left, donor_id)
        right_reactive = self._drop_one_chebi(right, donor_product_id)
        left_reactive = [
            participant for participant in left_reactive if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_reactive = [
            participant for participant in right_reactive if participant.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS}
        ]

        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive nitrogenous acceptor remains after removing donor phosphates",
            )

        if not any(self._has_nitrogen(participant) for participant in left_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No nitrogenous acceptor substrate detected",
            )

        left_total_p = sum(self._phosphorus_count(participant) for participant in left_reactive)
        right_total_p = sum(self._phosphorus_count(participant) for participant in right_reactive)
        if right_total_p != left_total_p + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Acceptor-side phosphorus count does not increase by one",
            )

        left_pn_bonds = self._pattern_count(left_reactive, self.PHOSPHORAMIDATE_PATTERN)
        right_pn_bonds = self._pattern_count(right_reactive, self.PHOSPHORAMIDATE_PATTERN)
        if right_pn_bonds <= left_pn_bonds:
            return ClassificationResult(
                is_member=False,
                explanation="Reaction does not form a new phosphorus-nitrogen bond",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Phosphotransferase with nitrogenous-group acceptor: a phosphorus donor "
                "phosphorylates a nitrogen center to form a phosphoramidate"
            ),
        )

    @classmethod
    def _find_donor_pair(
        cls,
        left: list[Participant],
        right: list[Participant],
    ) -> tuple[str, str] | None:
        left_ids = {participant.chebi_id for participant in left if participant.chebi_id}
        right_ids = {participant.chebi_id for participant in right if participant.chebi_id}
        for donor_id, donor_product_id in cls.DONOR_PRODUCT_PAIRS.items():
            if donor_id in left_ids and donor_product_id in right_ids:
                return donor_id, donor_product_id
        return None

    @staticmethod
    def _drop_one_chebi(participants: list[Participant], chebi_id: str) -> list[Participant]:
        """Remove one instance of a donor or spent donor."""
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
    def _has_nitrogen(participant: Participant) -> bool:
        """Check whether a participant has at least one nitrogen atom."""
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        """Count phosphorus atoms in a participant."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    @classmethod
    def _pattern_count(
        cls,
        participants: list[Participant],
        pattern: Chem.Mol | None,
    ) -> int:
        """Count occurrences of a SMARTS pattern across participants."""
        if pattern is None:
            return 0
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total
