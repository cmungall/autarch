"""Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or electrons are transferred from one donor, and one oxygen atom is incorporated into a donor.

This EC 1.13.12 branch covers internal monooxygenases and internal mixed
function oxidases. These reactions use dioxygen directly, do not require an
external reductant such as NAD(P)H or flavin, and typically release carbon
dioxide, water, or light as the second oxygen-derived outcome.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_DIPHOSPHATE,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_FMN,
    CHEBI_FMNH2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
CHEBI_SUCCINATE = "CHEBI:30031"
CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"
CHEBI_FE2 = "CHEBI:29033"
CHEBI_FE3 = "CHEBI:29034"


class OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfOneAtomOfOxygenInternalMonooxygenasesOrInternalMixedFunctionOxidases(Oxidoreductase):
    """oxidoreductase acting on single donors with incorporation of molecular oxygen, incorporation of one atom of oxygen (internal monooxygenases or internal mixed function oxidases)

    Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or
    electrons are transferred from one donor, and one oxygen atom is
    incorporated into a donor.

    The defining signature is a single substantive donor plus dioxygen, with no
    external nicotinamide, flavin, ferredoxin, or 2-oxoglutarate cochemistry.
    ATP may be used for luciferin activation, but the oxygenation itself remains
    internally driven.
    """

    GO_ID = "GO:0016703"
    EC_NUMBER_PREFIX = "1.13.12.-"

    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
        CHEBI_H2O2,
        CHEBI_FERREDOXIN_OXIDIZED,
        CHEBI_FERREDOXIN_REDUCED,
        CHEBI_FE2,
        CHEBI_FE3,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {
        CHEBI_CO2,
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_AMP,
        CHEBI_ADP,
        CHEBI_PHOSPHATE,
        CHEBI_DIPHOSPHATE,
    }
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for internal monooxygenase chemistry."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not internal monooxygenase chemistry",
            )

        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
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
                explanation="External nicotinamide, flavin, ferredoxin, peroxide, or 2-oxoglutarate cochemistry indicates another oxygenase branch",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen cosubstrate detected",
            )

        if any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Hydrogen peroxide production indicates an oxidase rather than internal monooxygenase chemistry",
            )

        if not any(participant.chebi_id == CHEBI_CO2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Internal monooxygenases in this branch are treated as decarboxylative oxygenations and should release carbon dioxide",
            )

        if CHEBI_ATP in chebis:
            if not any(participant.chebi_id == CHEBI_AMP for participant in right_participants):
                return ClassificationResult(
                    is_member=False,
                    explanation="ATP-dependent internal monooxygenases should convert ATP to AMP",
                )
            if not any(participant.chebi_id == CHEBI_DIPHOSPHATE for participant in right_participants):
                return ClassificationResult(
                    is_member=False,
                    explanation="ATP-dependent internal monooxygenases should release diphosphate",
                )
        if self._count_chebi(right_participants, CHEBI_H2O) > 1:
            return ClassificationResult(
                is_member=False,
                explanation="Formation of multiple water molecules indicates a different oxygenase branch",
            )

        if self._has_two_oxoglutarate_succinate_branch(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="2-oxoglutarate/succinate cochemistry indicates a different oxygenase branch",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS and self._has_carbon(participant)
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS and self._has_carbon(participant)
        ]
        if len(left_reactive) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="This branch acts on a single substantive donor after removing dioxygen and activation cofactors",
            )
        if not right_reactive:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive organic product remains after removing gas and activation byproducts",
            )

        if CHEBI_ATP not in chebis and not any(
            participant.chebi_id == CHEBI_H2O for participant in right_participants
        ):
            if not self._has_amide_gain(left_reactive, right_reactive):
                return ClassificationResult(
                    is_member=False,
                    explanation="Non-ATP internal monooxygenases without water release should show amide formation in the principal product",
                )

        left_carbons = sum(self._carbon_count(participant) for participant in left_reactive)
        right_carbons = sum(self._carbon_count(participant) for participant in right_reactive)
        if not self._is_two_oxoglutarate_self_oxidation(left_reactive, right_participants):
            if left_carbons - right_carbons not in {0, 1}:
                return ClassificationResult(
                    is_member=False,
                    explanation="The principal product scaffold does not match a single-donor oxygenation or oxidative decarboxylation",
                )

        return ClassificationResult(
            is_member=True,
            explanation="Internal monooxygenase/internal mixed function oxidase: a single donor reacts with dioxygen without an external reductant",
        )

    @classmethod
    def _has_amide_gain(
        cls,
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        return cls._substructure_count(right_reactive, cls.AMIDE_PATTERN) > cls._substructure_count(
            left_reactive,
            cls.AMIDE_PATTERN,
        )

    @staticmethod
    def _is_two_oxoglutarate_self_oxidation(
        left_reactive: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        if len(left_reactive) != 1 or left_reactive[0].chebi_id != CHEBI_2_OXOGLUTARATE:
            return False
        return sum(
            max(participant.count, 1)
            for participant in right_participants
            if participant.chebi_id == CHEBI_CO2
        ) >= 3

    @staticmethod
    def _substructure_count(
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
    def _count_chebi(participants: list[Participant], chebi_id: str) -> int:
        return sum(
            max(participant.count, 1)
            for participant in participants
            if participant.chebi_id == chebi_id
        )

    @staticmethod
    def _has_two_oxoglutarate_succinate_branch(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        return any(
            participant.chebi_id == CHEBI_2_OXOGLUTARATE
            for participant in left_participants
        ) and any(
            participant.chebi_id == CHEBI_SUCCINATE
            for participant in right_participants
        )

    @staticmethod
    def _has_carbon(participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    @staticmethod
    def _carbon_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) * max(participant.count, 1)
