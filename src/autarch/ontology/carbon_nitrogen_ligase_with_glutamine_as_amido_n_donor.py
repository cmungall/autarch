"""Catalysis of the transfer of the amide nitrogen of glutamine to a substrate. Usually composed of two subunits or domains, one that first hydrolyzes glutamine, and then transfers the resulting ammonia to the second subunit (or domain), where it acts as a source of nitrogen.

This EC 6.3.5 branch captures ATP- or NTP-coupled glutamine amidotransferase
ligation chemistry. The common mechanistic signature is glutamine hydrolysis to
L-glutamate coupled to productive nitrogen transfer onto another substrate.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.ligase import Ligase

CHEBI_GLUTAMINE = "CHEBI:58359"
CHEBI_GLUTAMATE = "CHEBI:29985"
CHEBI_POLYPHOSPHATE = "CHEBI:16838"


class CarbonNitrogenLigaseWithGlutamineAsAmidoNDonor(Ligase):
    """carbon-nitrogen ligase with glutamine as amido-N-donor

    Catalysis of the transfer of the amide nitrogen of glutamine to a
    substrate. Usually composed of two subunits or domains, one that first
    hydrolyzes glutamine, and then transfers the resulting ammonia to the
    second subunit (or domain), where it acts as a source of nitrogen.

    The class is defined by glutamine-to-glutamate conversion coupled to
    triphosphate consumption and a substantive change in another substrate that
    receives the transferred amido nitrogen.
    """

    GO_ID = "GO:0016884"
    EC_NUMBER_PREFIX = "6.3.5.-"

    NTP_REACTANTS = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}
    SPENT_NUCLEOTIDES = {
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_GDP,
        CHEBI_UDP,
        CHEBI_UMP,
        CHEBI_CDP,
        CHEBI_CMP,
    }
    EXCLUDED_CHEBIS = {
        CHEBI_O2,
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
    }
    LEFT_SPECTATOR_CHEBIS = NTP_REACTANTS | {CHEBI_GLUTAMINE, CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = SPENT_NUCLEOTIDES | {
        CHEBI_GLUTAMATE,
        CHEBI_PHOSPHATE,
        CHEBI_DIPHOSPHATE,
        CHEBI_POLYPHOSPHATE,
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    PHOSPHORUS_NITROGEN_PATTERN = Chem.MolFromSmarts("[PX4](=[OX1])([O])([O])[NX3]")
    POLYMER_TYPES = {PolymerType.PROTEIN, PolymerType.PEPTIDE, PolymerType.TRNA, PolymerType.RNA}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for glutamine-dependent carbon-nitrogen ligase chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not glutamine-dependent ligase chemistry",
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
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External redox or oxygen chemistry indicates another ligase branch",
            )

        if not any(participant.chebi_id == CHEBI_GLUTAMINE for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No glutamine amide donor detected",
            )

        if not any(participant.chebi_id == CHEBI_GLUTAMATE for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No glutamate coproduct detected from glutamine hydrolysis",
            )

        if not any(participant.chebi_id == CHEBI_H2O for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Glutamine-dependent amidotransfer ligases consume water during the glutaminase half-reaction",
            )

        if not self._has_triphosphate_coupling(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No ATP-or-similar triphosphate coupling detected",
            )

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
        if not left_reactive or not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive acceptor/product pair remains after removing glutamine, glutamate, and nucleotide spectators",
            )

        if not self._has_productive_n_transfer(left_reactive, right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="No productive nitrogen-transfer signature remains after glutamine hydrolysis",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Carbon-nitrogen ligase with glutamine as amido-N donor: glutamine hydrolysis is coupled to triphosphate-driven nitrogen transfer",
        )

    @classmethod
    def _has_triphosphate_coupling(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        has_ntp = any(participant.chebi_id in cls.NTP_REACTANTS for participant in left_participants)
        has_spent_nucleotide = any(
            participant.chebi_id in cls.SPENT_NUCLEOTIDES for participant in right_participants
        )
        has_phosphate_product = any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, CHEBI_POLYPHOSPHATE}
            for participant in right_participants
        )
        return has_ntp and has_spent_nucleotide and has_phosphate_product

    @classmethod
    def _has_productive_n_transfer(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        left_amides = cls._substructure_count(left_reactive, cls.AMIDE_PATTERN)
        right_amides = cls._substructure_count(right_reactive, cls.AMIDE_PATTERN)
        if right_amides > left_amides:
            return True

        left_pn = cls._substructure_count(left_reactive, cls.PHOSPHORUS_NITROGEN_PATTERN)
        right_pn = cls._substructure_count(right_reactive, cls.PHOSPHORUS_NITROGEN_PATTERN)
        if right_pn > left_pn:
            return True

        if cls._has_polymer_conversion(left_reactive, right_reactive):
            return True

        left_ids = {cls._participant_key(participant) for participant in left_reactive}
        right_ids = {cls._participant_key(participant) for participant in right_reactive}
        return left_ids != right_ids

    @classmethod
    def _has_polymer_conversion(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        left_polymers = [p for p in left_reactive if p.polymer_type in cls.POLYMER_TYPES]
        right_polymers = [p for p in right_reactive if p.polymer_type in cls.POLYMER_TYPES]
        return any(
            left_polymer.polymer_type == right_polymer.polymer_type
            and not left_polymer.is_same_molecule(right_polymer)
            for left_polymer in left_polymers
            for right_polymer in right_polymers
        )

    @classmethod
    def _substructure_count(
        cls,
        participants: list[Participant],
        pattern: Chem.Mol | None,
    ) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None or pattern is None:
                continue
            total += len(mol.GetSubstructMatches(pattern)) * max(participant.count, 1)
        return total

    @staticmethod
    def _participant_key(participant: Participant) -> tuple[str | None, str | None, PolymerType | None]:
        return participant.chebi_id, participant.smiles, participant.polymer_type
