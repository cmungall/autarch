"""oxidoreductase acting on the aldehyde or oxo group of donors, oxygen as
acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde or
ketone (oxo) group acts as a hydrogen or electron donor and reduces oxygen.
"""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety, is_aldehyde, is_carboxylic_acid, is_ketone
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_FAD,
    CHEBI_FADH2,
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


class OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase acting on the aldehyde or oxo group of donors, oxygen as
    acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde
    or ketone (oxo) group acts as a hydrogen or electron donor and reduces
    oxygen.
    """

    GO_ID = "GO:0016623"
    EC_NUMBER_PREFIX = "1.2.3.-"

    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        "CHEBI:58210",  # FMN
        "CHEBI:58307",  # FMNH2
        "CHEBI:33722",  # oxidized [4Fe-4S]
        "CHEBI:33723",  # reduced [4Fe-4S]
        "CHEBI:33737",  # oxidized ferredoxin
        "CHEBI:33738",  # reduced ferredoxin
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O2, CHEBI_H2O, CHEBI_H_PLUS}
    AMMONIA_CHEBIS = {CHEBI_NH3, CHEBI_NH4}
    FORMATE_PATTERN = Chem.MolFromSmarts("[CH1](=O)[O;H1,X1-]")
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=O)O[P](=O)")
    CARBOXYL_PATTERN = Chem.MolFromSmarts("[CX3](=O)[O;H1,X1-]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for aldehyde/oxo oxidation with oxygen as the terminal acceptor."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not aldehyde/oxo oxidoreductase chemistry",
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
        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing molecular oxygen acceptor on the substrate side",
            )
        if not any(participant.chebi_id == CHEBI_H2O2 for participant in right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Missing hydrogen peroxide product expected for oxygen-accepting aldehyde/oxo redox",
            )
        chebis = {
            participant.chebi_id
            for participant in left_participants + right_participants
            if participant.chebi_id
        }
        if chebis & self.EXCLUDED_CHEBIS:
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide, flavin, or iron-sulfur cofactors indicate a different oxidoreductase branch",
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
                explanation="No substantive donor/product pair remained after removing oxygen and solvent spectators",
            )

        if any(participant.chebi_id in self.AMMONIA_CHEBIS for participant in right_reactive):
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia-releasing oxygen redox belongs to CH-N donor oxidoreductase branches",
            )

        donor_candidates = [
            participant for participant in left_reactive if self._is_oxo_donor(participant)
        ]
        if not donor_candidates:
            return ClassificationResult(
                is_member=False,
                explanation="No aldehyde or oxo donor detected on the substrate side",
            )

        if any(self._is_oxidized_product(participant) for participant in right_reactive):
            return ClassificationResult(
                is_member=True,
                explanation="Aldehyde/oxo oxidoreductase: oxygen-dependent oxidation to a carboxylate, thioester, or acyl phosphate",
            )

        co2_products = sum(
            max(participant.count, 1)
            for participant in right_reactive
            if participant.chebi_id == CHEBI_CO2
        )
        if co2_products >= 1 and any(
            self._supports_decarboxylative_branch(participant)
            for participant in donor_candidates
        ):
            return ClassificationResult(
                is_member=True,
                explanation="Aldehyde/oxo oxidoreductase: oxygen-dependent oxidative decarboxylation of an oxo acid donor",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No oxidized product matched aldehyde/oxo donor oxygen-acceptor chemistry",
        )

    @classmethod
    def _is_oxo_donor(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        if is_aldehyde(participant.smiles, participant.chebi_id) or cls._is_formate_like(participant):
            return True
        if is_ketone(participant.smiles, participant.chebi_id):
            return True
        return cls._is_oxalate_like(participant)

    @staticmethod
    def _is_oxidized_product(participant: Participant) -> bool:
        return bool(
            (participant.smiles and is_carboxylic_acid(participant.smiles, participant.chebi_id))
            or participant.is_thioester()
            or OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsOxygenAsAcceptor._is_acyl_phosphate_product(participant)
        )

    @classmethod
    def _supports_decarboxylative_branch(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        return (
            is_ketone(participant.smiles, participant.chebi_id)
            and participant.has_moiety(Moiety.CARBOXYL)
        ) or cls._is_oxalate_like(participant)

    @classmethod
    def _is_formate_like(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        return mol is not None and cls.FORMATE_PATTERN is not None and mol.HasSubstructMatch(cls.FORMATE_PATTERN)

    @classmethod
    def _is_acyl_phosphate_product(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        return mol is not None and cls.ACYL_PHOSPHATE_PATTERN is not None and mol.HasSubstructMatch(cls.ACYL_PHOSPHATE_PATTERN)

    @classmethod
    def _carboxyl_group_count(cls, participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None or cls.CARBOXYL_PATTERN is None:
            return 0
        return len(mol.GetSubstructMatches(cls.CARBOXYL_PATTERN))

    @classmethod
    def _is_oxalate_like(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        hetero_n_or_s = any(atom.GetAtomicNum() in {7, 16} for atom in mol.GetAtoms())
        return cls._carboxyl_group_count(participant) >= 2 and carbon_count <= 2 and not hetero_n_or_s
