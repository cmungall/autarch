"""sulfurtransferase activity.

Catalysis of the transfer of sulfur atoms from one compound (donor) to another
(acceptor).
"""

from collections import Counter

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.transferase import Transferase

CHEBI_PAPS = "CHEBI:58339"
CHEBI_PAP = "CHEBI:58343"


class Sulfurtransferase(Transferase):
    """sulfurtransferase activity.

    Catalysis of the transfer of sulfur atoms from one compound (donor) to
    another (acceptor).
    """

    GO_ID = "GO:0016783"
    EC_NUMBER_PREFIX = "2.8.1.-"
    LEFT_SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_ATP,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_NADH,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_O2,
        CHEBI_H2O2,
    }
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_DIPHOSPHATE,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_NADH,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_H2O2,
    }
    REDOX_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_NADH,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_O2,
        CHEBI_H2O2,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for sulfur transfer without nicotinamide/flavin redox chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sulfur transfer",
            )

        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        if chebis & self.REDOX_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="Nicotinamide, flavin, or oxygen cofactors indicate a different sulfur chemistry branch",
            )
        if CHEBI_SAM in chebis or CHEBI_SAH in chebis:
            return ClassificationResult(
                is_member=False,
                explanation="SAM/SAH chemistry indicates sulfur-containing group transfer rather than sulfur atom transfer",
            )
        if CHEBI_PAPS in chebis or CHEBI_PAP in chebis:
            return ClassificationResult(
                is_member=False,
                explanation="PAPS/PAP chemistry indicates sulfotransferase rather than sulfurtransferase activity",
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
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS and participant.get_mol() is not None
        ]
        if len(left_reactive) < 2 or len(right_reactive) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Sulfur transfer requires a donor, an acceptor, and chemically resolved products",
            )

        donor_pairs = self._find_desulfurized_pairs(left_reactive, right_reactive)
        if not donor_pairs:
            return ClassificationResult(
                is_member=False,
                explanation="No sulfur donor loses sulfur while retaining its non-sulfur scaffold",
            )

        if not self._has_transfer_signature(left_reactive, right_reactive, donor_pairs):
            return ClassificationResult(
                is_member=False,
                explanation="No complementary sulfur transfer signature detected in the products",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sulfurtransferase: sulfur donor loses sulfur while another product retains the transferred sulfur",
        )

    def _find_desulfurized_pairs(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> list[tuple[Participant, Participant]]:
        donor_pairs: list[tuple[Participant, Participant]] = []
        for substrate in left_reactive:
            substrate_s = self._sulfur_count(substrate)
            if substrate_s == 0:
                continue
            substrate_backbone = self._sulfur_stripped_scaffold(substrate)
            substrate_carbons = self._carbon_count(substrate)
            for product in right_reactive:
                product_s = self._sulfur_count(product)
                if product_s >= substrate_s:
                    continue
                product_backbone = self._sulfur_stripped_scaffold(product)
                if substrate_carbons == 0 and self._carbon_count(product) == 0:
                    non_s_counts = self._non_sulfur_heavy_atom_counts(substrate)
                    if non_s_counts and non_s_counts == self._non_sulfur_heavy_atom_counts(product):
                        donor_pairs.append((substrate, product))
                elif substrate_backbone and substrate_backbone == product_backbone:
                    donor_pairs.append((substrate, product))
        return donor_pairs

    def _has_transfer_signature(
        self,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
        donor_pairs: list[tuple[Participant, Participant]],
    ) -> bool:
        matched_left = {id(left) for left, _ in donor_pairs}
        matched_right = {id(right) for _, right in donor_pairs}

        for substrate in left_reactive:
            if id(substrate) in matched_left:
                continue
            substrate_s = self._sulfur_count(substrate)
            substrate_backbone = self._sulfur_stripped_scaffold(substrate)
            if not substrate_backbone:
                continue
            for product in right_reactive:
                if id(product) in matched_right:
                    continue
                if self._sulfur_count(product) <= substrate_s:
                    continue
                if substrate_backbone == self._sulfur_stripped_scaffold(product):
                    return True

        unmatched_sulfur_products = [
            participant
            for participant in right_reactive
            if id(participant) not in matched_right and self._sulfur_count(participant) > 0
        ]
        return bool(unmatched_sulfur_products)

    @staticmethod
    def _sulfur_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    @staticmethod
    def _non_sulfur_heavy_atom_counts(participant: Participant) -> Counter[int]:
        mol = participant.get_mol()
        if mol is None:
            return Counter()
        return Counter(
            atom.GetAtomicNum()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() not in {1, 16}
        )

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    @staticmethod
    def _sulfur_stripped_scaffold(participant: Participant) -> str:
        mol = participant.get_mol()
        if mol is None:
            return ""
        sulfur = Chem.MolFromSmarts("[#16]")
        if sulfur is None:
            return ""
        stripped = Chem.DeleteSubstructs(Chem.Mol(mol), sulfur)
        if stripped.GetNumAtoms() == 0:
            return ""
        Chem.RemoveStereochemistry(stripped)
        return Chem.MolToSmiles(stripped, canonical=True, isomericSmiles=False)
