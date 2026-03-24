"""hexosyltransferase.

Catalysis of the transfer of a hexosyl group from one compound (donor) to
another (acceptor).
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
    CHEBI_CTP,
)
from autarch.ontology.glycosyltransferase import Glycosyltransferase

CHEBI_GMP = "CHEBI:58115"


class Hexosyltransferase(Glycosyltransferase):
    """hexosyltransferase.

    Catalysis of the transfer of a hexosyl group from one compound (donor) to
    another (acceptor).
    """

    GO_ID = "GO:0016758"
    EC_NUMBER_PREFIX = "2.4.1.-"
    NUCLEOTIDE_HEXOSE_DONOR_IDS = {
        "CHEBI:18066",  # UDP-D-glucose
        "CHEBI:18307",  # UDP-D-galactose
        "CHEBI:58885",  # UDP-alpha-D-glucose(2-)
        "CHEBI:66914",  # UDP-alpha-D-galactose(2-)
        "CHEBI:57705",  # UDP-N-acetyl-alpha-D-glucosamine
        "CHEBI:58052",  # UDP-alpha-D-glucuronate
        "CHEBI:67138",  # UDP-N-acetyl-alpha-D-galactosamine
        "CHEBI:68623",  # UDP-N-acetyl-alpha-D-mannosamine
        "CHEBI:83836",  # UDP-rhamnose
        "CHEBI:57527",  # GDP-mannose
        "CHEBI:57273",  # GDP-fucose
        "CHEBI:62230",  # GDP-glucose
        "CHEBI:57498",  # ADP-alpha-D-glucose
        "CHEBI:61506",  # CDP-linked hexose
        "CHEBI:76533",  # NDP-alpha-D-glucose
        "CHEBI:57930",  # dTDP-linked deoxyhexose
        "CHEBI:70731",  # UDP-ManNAcA
        "CHEBI:70784",  # CDP-abequose
    }
    FREE_NDP_PRODUCTS = {CHEBI_UDP, CHEBI_GDP, CHEBI_ADP, CHEBI_CDP}
    FREE_NMP_PRODUCTS = {CHEBI_UMP, CHEBI_CMP, CHEBI_GMP}
    NTP_REACTANTS = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for EC 2.4.1-style hexosyl transfer chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not hexosyl transfer",
            )

        forward_result = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward_result.is_member:
            return forward_result

        reverse_result = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse_result.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse_result.explanation} (reverse reaction orientation)",
            )

        return forward_result

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        if self._is_hexose_donor_synthesis(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide-sugar donor synthesis, not hexosyl transfer",
            )

        if self._is_hexose_donor_interconversion(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide-sugar interconversion, not transfer to an external acceptor",
            )

        if self._is_phosphoglycosyl_transferase(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Phosphoglycosyl transferase chemistry transfers hexose-phosphate, not hexose",
            )

        if self._is_nucleotide_hexosyl_transfer(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Hexosyltransferase: nucleotide-hexose donor transfers a hexosyl group to an acceptor",
            )

        if self._is_hexose_phosphorylase(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Hexosyltransferase: phosphate accepts a transferred hexosyl group",
            )

        if self._is_carbohydrate_transglycosylation(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Hexosyltransferase: carbohydrate donor transfers a hexosyl group to another acceptor",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No hexosyl-transfer pattern detected",
        )

    @classmethod
    def _is_hexose_donor_synthesis(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect nucleotide-sugar biosynthesis rather than glycosyl transfer."""
        has_ntp = any(participant.chebi_id in cls.NTP_REACTANTS for participant in left_participants)
        has_hexose_phosphate = any(cls._is_hexose_phosphate(participant) for participant in left_participants)
        has_nucleotide_hexose = any(
            participant.chebi_id in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS for participant in right_participants
        )
        has_phosphate_byproduct = any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE}
            for participant in right_participants
        )
        return has_ntp and has_hexose_phosphate and has_nucleotide_hexose and has_phosphate_byproduct

    @classmethod
    def _is_phosphoglycosyl_transferase(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect phosphoglycosyl transferases that release NMP rather than NDP."""
        has_hexose_donor = any(
            participant.chebi_id in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS for participant in left_participants
        )
        has_nmp_product = any(
            participant.chebi_id in cls.FREE_NMP_PRODUCTS for participant in right_participants
        )
        return has_hexose_donor and has_nmp_product

    @classmethod
    def _is_nucleotide_hexosyl_transfer(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect direct transfer from a nucleotide-hexose donor."""
        has_hexose_donor = any(
            participant.chebi_id in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS for participant in left_participants
        )
        has_ndp_product = any(
            participant.chebi_id in cls.FREE_NDP_PRODUCTS for participant in right_participants
        )
        has_acceptor = any(
            participant.chebi_id not in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS
            and participant.chebi_id not in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE}
            for participant in left_participants
        )
        has_hexosylated_product = any(
            participant.chebi_id not in cls.FREE_NDP_PRODUCTS and cls._is_hexose_related(participant)
            for participant in right_participants
        )
        return has_hexose_donor and has_ndp_product and has_acceptor and has_hexosylated_product

    @classmethod
    def _is_hexose_phosphorylase(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect phosphorolysis/synthesis within the EC 2.4.1 branch."""
        has_phosphate = any(
            participant.chebi_id == CHEBI_PHOSPHATE for participant in left_participants
        )
        if not has_phosphate:
            return False
        if any(participant.chebi_id in cls.NTP_REACTANTS for participant in left_participants + right_participants):
            return False
        has_hexose_donor = any(cls._is_hexose_related(participant) for participant in left_participants)
        has_hexose_phosphate = any(
            cls._is_hexose_phosphate(participant) for participant in right_participants
        )
        has_hexose_product = any(
            cls._is_hexose_related(participant) for participant in right_participants
        )
        return has_hexose_donor and has_hexose_phosphate and has_hexose_product

    @classmethod
    def _is_carbohydrate_transglycosylation(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect transfer between carbohydrate or glycoside donors/acceptors."""
        if any(
            participant.chebi_id in cls.NTP_REACTANTS | cls.FREE_NDP_PRODUCTS | cls.FREE_NMP_PRODUCTS
            for participant in left_participants + right_participants
        ):
            return False
        if any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE}
            for participant in left_participants + right_participants
        ):
            return False
        if len(left_participants) == 1 and len(right_participants) == 1:
            return False

        left_hexose_like = sum(1 for participant in left_participants if cls._is_hexose_related(participant))
        right_hexose_like = sum(1 for participant in right_participants if cls._is_hexose_related(participant))
        return left_hexose_like >= 2 and right_hexose_like >= 2

    @classmethod
    def _is_hexose_donor_interconversion(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Detect nucleotide-sugar synthases and interconversions within the donor pool."""
        allowed_chebis = (
            cls.NUCLEOTIDE_HEXOSE_DONOR_IDS
            | cls.FREE_NDP_PRODUCTS
            | cls.FREE_NMP_PRODUCTS
            | cls.NTP_REACTANTS
            | {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE}
        )
        all_participants = left_participants + right_participants
        if not all(
            participant.chebi_id in allowed_chebis or cls._is_hexose_phosphate(participant)
            for participant in all_participants
        ):
            return False
        has_donor = any(
            participant.chebi_id in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS for participant in all_participants
        )
        has_hexose_phosphate = any(cls._is_hexose_phosphate(participant) for participant in all_participants)
        return has_donor and has_hexose_phosphate

    @classmethod
    def _is_hexose_related(cls, participant: Participant) -> bool:
        """Check for free or conjugated hexose/hexuronate-like participation."""
        if participant.chebi_id in cls.NUCLEOTIDE_HEXOSE_DONOR_IDS:
            return True
        if cls._has_hexose_like_ring(participant):
            return True

        return cls._is_hexose_phosphate(participant)

    @classmethod
    def _is_hexose_phosphate(cls, participant: Participant) -> bool:
        """Check for a phosphorylated hexose-like species."""
        counts = cls._element_counts(participant)
        if counts["P"] == 0 or counts["S"] != 0:
            return False
        return cls._has_hexose_like_ring(participant)

    @staticmethod
    def _element_counts(participant: Participant) -> dict[str, int]:
        """Count elements in a participant."""
        mol = participant.get_mol()
        counts = {"C": 0, "O": 0, "N": 0, "P": 0, "S": 0}
        if mol is None:
            return counts
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in counts:
                counts[symbol] += 1
        return counts

    @staticmethod
    def _has_hexose_like_ring(participant: Participant) -> bool:
        """Check for a six-membered oxygen-containing sugar ring."""
        mol = participant.get_mol()
        if mol is None:
            return False
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) != 6:
                continue
            ring_atoms = [mol.GetAtomWithIdx(index) for index in ring]
            oxygen_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            carbon_in_ring = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if oxygen_in_ring != 1 or carbon_in_ring < 5:
                continue
            hetero_substituents = 0
            ring_atom_ids = set(ring)
            for atom_index in ring:
                atom = mol.GetAtomWithIdx(atom_index)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring_atom_ids:
                        continue
                    if neighbor.GetAtomicNum() in {7, 8, 15}:
                        hetero_substituents += 1
                        break
            if hetero_substituents >= 2:
                return True
        return False
