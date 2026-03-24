"""phosphotransferase activity, carboxyl group as acceptor.

Catalysis of the transfer of a phosphorus-containing group from one compound
(donor) to a carboxyl group (acceptor).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
)
from autarch.moiety import is_carboxylic_acid
from autarch.ontology.transferase import Transferase


class PhosphotransferaseCarboxylGroupAsAcceptor(Transferase):
    """phosphotransferase activity, carboxyl group as acceptor.

    Catalysis of the transfer of a phosphorus-containing group from one
    compound (donor) to a carboxyl group (acceptor).
    """

    GO_ID = "GO:0016774"
    EC_NUMBER_PREFIX = "2.7.2.-"
    DONOR_PRODUCT_PAIRS = {
        CHEBI_ATP: CHEBI_ADP,
        CHEBI_GTP: CHEBI_GDP,
    }
    LEFT_SPECTATOR_CHEBIS = set(DONOR_PRODUCT_PAIRS)
    RIGHT_SPECTATOR_CHEBIS = set(DONOR_PRODUCT_PAIRS.values()) | {CHEBI_H2O, CHEBI_H_PLUS}
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])O[P]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for phosphoryl transfer onto a carboxyl group."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not phosphotransfer to a carboxyl group",
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
                explanation="No substantive acceptor substrate remains after removing donor phosphates",
            )

        if donor_id in {CHEBI_ATP, CHEBI_GTP} and any(
            participant.chebi_id in {CHEBI_PHOSPHATE, "CHEBI:16838"} for participant in right_reactive
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Free phosphate release indicates coupled ligase chemistry rather than direct carboxyl-group phosphotransfer",
            )

        left_acids = [participant for participant in left_reactive if self._is_carboxylate_like(participant)]
        if not left_acids:
            return ClassificationResult(
                is_member=False,
                explanation="No carboxyl-group acceptor substrate detected",
            )

        left_total_p = sum(self._phosphorus_count(participant) for participant in left_reactive)
        right_total_p = sum(self._phosphorus_count(participant) for participant in right_reactive)
        if right_total_p != left_total_p + 1:
            return ClassificationResult(
                is_member=False,
                explanation="Acceptor-side phosphorus count does not increase by one",
            )

        left_acyl_phosphates = self._pattern_count(left_reactive, self.ACYL_PHOSPHATE_PATTERN)
        right_acyl_phosphates = self._pattern_count(right_reactive, self.ACYL_PHOSPHATE_PATTERN)
        if right_acyl_phosphates <= left_acyl_phosphates:
            return ClassificationResult(
                is_member=False,
                explanation="Reaction does not form a new acyl-phosphate linkage",
            )

        return ClassificationResult(
            is_member=True,
            explanation=(
                "Phosphotransferase with carboxyl-group acceptor: a phosphorus donor "
                "phosphorylates a carboxylate to form an acyl phosphate"
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
        if CHEBI_DIPHOSPHATE in left_ids and ({CHEBI_PHOSPHATE, "CHEBI:16838"} & right_ids):
            return CHEBI_DIPHOSPHATE, "CHEBI:16838" if "CHEBI:16838" in right_ids else CHEBI_PHOSPHATE
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

    @staticmethod
    def _is_carboxylate_like(participant: Participant) -> bool:
        """Check whether a participant contains a carboxyl group."""
        return bool(participant.smiles and is_carboxylic_acid(participant.smiles, participant.chebi_id))
