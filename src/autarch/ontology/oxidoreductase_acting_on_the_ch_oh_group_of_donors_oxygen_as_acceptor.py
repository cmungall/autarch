"""oxidoreductase acting on the CH-OH group of donors oxygen as acceptor.

Catalysis of an oxidation-reduction reaction in which the donor is a CH-OH group
and molecular oxygen is the acceptor.
"""

import json
from pathlib import Path

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety, has_moiety, is_aldehyde, is_carboxylic_acid, is_ketone
from autarch.molecules import (
    CHEBI_COA,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase


class OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on the CH-OH group of donors oxygen as acceptor.

    Catalysis of an oxidation-reduction reaction in which the donor is a CH-OH
    group and molecular oxygen is the acceptor.
    """

    GO_ID = "GO:0016899"  # oxidoreductase activity, acting on CH-OH, oxygen as acceptor
    EC_NUMBER_PREFIX = "1.1.3.-"  # EC prefix for alcohol oxidases
    SOLUBLE_REDOX_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
    }
    SPECTATOR_CHEBIS = {
        CHEBI_COA,
        CHEBI_O2,
        CHEBI_H2O2,
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }
    AMMONIA_CHEBIS = {CHEBI_NH3, CHEBI_NH4}
    _CHEBI_SMILES_CACHE: dict[str, str] | None = None

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for O2-dependent oxidation of a hydroxyl-bearing donor."""
        has_oxygen = any(
            participant.chebi_id == CHEBI_O2 or self._effective_smiles(participant) == "O=O"
            for participant in reaction.left_participants
        )
        if not has_oxygen:
            return ClassificationResult(
                is_member=False,
                explanation="No oxygen (O2) reactant found",
            )

        has_peroxide = any(
            participant.chebi_id == CHEBI_H2O2
            or self._effective_smiles(participant) in {"OO", "[H]OO[H]"}
            for participant in reaction.right_participants
        )
        if not has_peroxide:
            return ClassificationResult(
                is_member=False,
                explanation="No hydrogen peroxide (H2O2) product found",
            )

        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        has_soluble_redox_cofactor = any(
            participant.chebi_id in self.SOLUBLE_REDOX_CHEBIS
            for participant in reaction.all_participants()
        )
        if has_soluble_redox_cofactor:
            return ClassificationResult(
                is_member=False,
                explanation="Uses soluble NAD(P) cofactor - not oxygen-as-acceptor EC 1.1.3 chemistry",
            )

        produces_ammonia = any(
            participant.chebi_id in self.AMMONIA_CHEBIS
            for participant in reaction.right_participants
        )
        if produces_ammonia:
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia-releasing oxidation - not CH-OH donor chemistry",
            )

        left_non_spectators = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.SPECTATOR_CHEBIS
        ]
        right_non_spectators = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.SPECTATOR_CHEBIS
        ]

        has_hydroxyl_donor = any(
            self._is_hydroxyl_donor(participant) for participant in left_non_spectators
        )
        if not has_hydroxyl_donor:
            return ClassificationResult(
                is_member=False,
                explanation="No CH-OH donor detected among non-cofactor substrates",
            )

        has_oxidized_product = any(
            self._is_oxidized_product(participant)
            for participant in right_non_spectators
        )
        has_oxo_donor_left = any(
            self._is_oxo_donor(participant) for participant in left_non_spectators
        )
        has_strongly_oxidized_product = any(
            self._is_strongly_oxidized_product(participant)
            for participant in right_non_spectators
        )
        if has_oxo_donor_left and has_strongly_oxidized_product:
            return ClassificationResult(
                is_member=False,
                explanation="Aldehyde or oxo donor oxidation branch - not CH-OH donor chemistry",
            )

        if has_oxidized_product:
            return ClassificationResult(
                is_member=True,
                explanation="Alcohol oxidase: CH-OH donor oxidized with oxygen acceptor and H2O2 formation",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No oxidized product consistent with CH-OH donor oxidation",
        )

    @staticmethod
    def _is_hydroxyl_donor(participant: Participant) -> bool:
        smiles = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor._effective_smiles(
            participant
        )
        return bool(smiles and has_moiety(smiles, Moiety.HYDROXYL, participant.chebi_id))

    @staticmethod
    def _is_oxidized_product(participant: Participant) -> bool:
        smiles = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor._effective_smiles(
            participant
        )
        if not smiles:
            return False
        return bool(
            is_aldehyde(smiles, participant.chebi_id)
            or is_ketone(smiles, participant.chebi_id)
            or is_carboxylic_acid(smiles, participant.chebi_id)
            or has_moiety(smiles, Moiety.ACYL, participant.chebi_id)
        )

    @staticmethod
    def _is_oxo_donor(participant: Participant) -> bool:
        smiles = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor._effective_smiles(
            participant
        )
        return bool(
            smiles
            and (
                is_aldehyde(smiles, participant.chebi_id)
                or is_ketone(smiles, participant.chebi_id)
            )
        )

    @staticmethod
    def _is_strongly_oxidized_product(participant: Participant) -> bool:
        smiles = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor._effective_smiles(
            participant
        )
        return bool(
            smiles
            and (
                is_carboxylic_acid(smiles, participant.chebi_id)
                or participant.is_thioester()
            )
        )

    @classmethod
    def _effective_smiles(cls, participant: Participant) -> str | None:
        if participant.smiles:
            return participant.smiles
        if not participant.chebi_id:
            return None
        if cls._CHEBI_SMILES_CACHE is None:
            cache_path = Path(__file__).resolve().parents[3] / "cache" / "chebi_smiles.json"
            if cache_path.exists():
                with open(cache_path) as stream:
                    cls._CHEBI_SMILES_CACHE = json.load(stream)
            else:
                cls._CHEBI_SMILES_CACHE = {}
        return cls._CHEBI_SMILES_CACHE.get(participant.chebi_id)
