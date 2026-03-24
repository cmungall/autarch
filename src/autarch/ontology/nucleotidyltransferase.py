"""nucleotidyltransferase.

Catalysis of the transfer of a nucleotidyl group from one compound (donor) to
another (acceptor).
"""

from autarch.datamodel import ClassificationResult, Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GMP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.reaction import ReactionClass

CHEBI_POLYPHOSPHATE = "CHEBI:16838"


class Nucleotidyltransferase(ReactionClass):
    """nucleotidyltransferase.

    Catalysis of the transfer of a nucleotidyl group from one compound (donor)
    to another (acceptor).
    """

    GO_ID = "GO:0016779"
    EC_NUMBER_PREFIX = "2.7.7.-"
    DONOR_IDS = {
        CHEBI_ATP,
        CHEBI_GTP,
        CHEBI_UTP,
        CHEBI_CTP,
        CHEBI_ADP,
        CHEBI_GDP,
        CHEBI_UDP,
        CHEBI_CDP,
    }
    FREE_NUCLEOTIDE_IDS = DONOR_IDS | {CHEBI_AMP, CHEBI_CMP, CHEBI_GMP, CHEBI_UMP}
    BYPRODUCT_IDS = {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, CHEBI_POLYPHOSPHATE}
    SPECTATOR_IDS = {CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for nucleotidyl transfer by mass-balanced donor/acceptor logic."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not nucleotidyl transfer",
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
        """Evaluate one transfer orientation."""
        if self._is_polymer_extension(left, right):
            return ClassificationResult(
                is_member=True,
                explanation="Nucleotidyltransferase: polymer growth by nucleotidyl transfer with phosphate byproduct release",
            )

        donors = [participant for participant in left if participant.chebi_id in self.DONOR_IDS]
        byproducts = [participant for participant in right if participant.chebi_id in self.BYPRODUCT_IDS]
        if not donors:
            return ClassificationResult(
                is_member=False,
                explanation="No nucleoside phosphate donor found",
            )
        if not byproducts:
            return ClassificationResult(
                is_member=False,
                explanation="No phosphate or diphosphate leaving group found",
            )
        if any(participant.chebi_id == CHEBI_AMP for participant in right):
            return ClassificationResult(
                is_member=False,
                explanation="AMP release indicates an ATP-dependent ligase or activation step, not nucleotidyl transfer",
            )

        left_reactive = [
            participant
            for participant in left
            if participant.chebi_id not in self.DONOR_IDS | self.SPECTATOR_IDS
        ]
        right_reactive = [
            participant
            for participant in right
            if participant.chebi_id not in self.BYPRODUCT_IDS | self.SPECTATOR_IDS
        ]
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive acceptor substrate or nucleotidylated product remains after removing donors and leaving groups",
            )

        for donor in donors:
            donor_carbons = self._carbon_count(donor)
            donor_phosphorus = self._phosphorus_count(donor)
            donor_nitrogen = self._nitrogen_count(donor)
            if donor_carbons == 0 or donor_phosphorus == 0 or donor_nitrogen == 0:
                continue

            for byproduct in byproducts:
                byproduct_phosphorus = self._phosphorus_count(byproduct)
                if byproduct_phosphorus == 0:
                    continue

                for acceptor in left_reactive:
                    if acceptor.is_polymer() or acceptor.get_mol() is None:
                        continue
                    acceptor_carbons = self._carbon_count(acceptor)
                    if acceptor_carbons == 0:
                        continue
                    acceptor_phosphorus = self._phosphorus_count(acceptor)
                    for product in right_reactive:
                        if product.is_polymer() or product.get_mol() is None:
                            continue
                        if product.chebi_id in self.FREE_NUCLEOTIDE_IDS:
                            continue
                        if self._carbon_count(product) != donor_carbons + acceptor_carbons:
                            continue
                        if self._phosphorus_count(product) != (
                            donor_phosphorus + acceptor_phosphorus - byproduct_phosphorus
                        ):
                            continue
                        if self._nitrogen_count(product) < donor_nitrogen:
                            continue
                        return ClassificationResult(
                            is_member=True,
                            explanation="Nucleotidyltransferase: a nucleoside phosphate donor transfers a nucleotidyl group to an acceptor with phosphate leaving-group release",
                        )

        return ClassificationResult(
            is_member=False,
            explanation="No nucleotidyl-transfer mass-balance signature detected",
        )

    def _is_polymer_extension(self, left: list[Participant], right: list[Participant]) -> bool:
        """Detect polymerase-like nucleotidyl transfer on polymer backbones."""
        left_polymers = [participant for participant in left if participant.is_polymer()]
        right_polymers = [participant for participant in right if participant.is_polymer()]
        if len(left_polymers) != 1 or len(right_polymers) != 1:
            return False

        left_polymer = left_polymers[0]
        right_polymer = right_polymers[0]
        if left_polymer.polymer_type != right_polymer.polymer_type:
            return False
        if left_polymer.polymer_type not in {PolymerType.RNA, PolymerType.DNA, PolymerType.POLYNUCLEOTIDE}:
            return False
        if left_polymer.polymer_index != "n" or right_polymer.polymer_index != "n+1":
            return False

        has_donor = any(participant.chebi_id in self.DONOR_IDS for participant in left)
        has_byproduct = any(participant.chebi_id in self.BYPRODUCT_IDS for participant in right)
        return has_donor and has_byproduct

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        """Count carbon atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        """Count nitrogen atoms."""
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _phosphorus_count(participant: Participant) -> int:
        """Count phosphorus atoms, using CHEBI when structures are absent for phosphate byproducts."""
        mol = participant.get_mol()
        if mol is not None:
            return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
        if participant.chebi_id == CHEBI_PHOSPHATE:
            return 1
        if participant.chebi_id in {CHEBI_DIPHOSPHATE, CHEBI_POLYPHOSPHATE}:
            return 2
        return 0
