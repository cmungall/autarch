"""oxidoreductase acting on the aldehyde or oxo group of donors NAD or NADP as
acceptor.

Catalysis of an oxidation-reduction reaction in which the donor is an aldehyde
or oxo group and NAD or NADP is the acceptor.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety, is_aldehyde, is_carboxylic_acid, is_ketone
from autarch.molecules import (
    CHEBI_CO2,
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


class OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on the aldehyde or oxo group of donors NAD or NADP
    as acceptor.

    Catalysis of an oxidation-reduction reaction in which the donor is an
    aldehyde or oxo group and NAD or NADP is the acceptor.
    """

    GO_ID = "GO:0016620"  # oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor
    EC_NUMBER_PREFIX = "1.2.-.-"
    ACCEPTOR_OXIDIZED_CHEBIS = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    ACCEPTOR_REDUCED_CHEBIS = {CHEBI_NADH, CHEBI_NADPH}
    SPECTATOR_CHEBIS = {
        CHEBI_H2O,
        CHEBI_H_PLUS,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
        CHEBI_NADH,
        CHEBI_NADPH,
    }
    AMMONIA_CHEBIS = {CHEBI_NH3, CHEBI_NH4}
    FORMATE_PATTERN = Chem.MolFromSmarts("[CH1](=O)[O;H1,X1-]")
    ACYL_PHOSPHATE_PATTERN = Chem.MolFromSmarts("[CX3](=O)O[P](=O)")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for EC 1.2 oxidation signatures without relying on labels."""
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        has_acceptor_left = any(
            participant.chebi_id in self.ACCEPTOR_OXIDIZED_CHEBIS
            for participant in reaction.left_participants
        )
        has_reduced_acceptor_right = any(
            participant.chebi_id in self.ACCEPTOR_REDUCED_CHEBIS
            for participant in reaction.right_participants
        )
        if not (has_acceptor_left and has_reduced_acceptor_right):
            return ClassificationResult(
                is_member=False,
                explanation="No NAD(P)+ acceptor -> NAD(P)H reduction signature",
            )

        has_oxygen_cochemistry = any(
            participant.chebi_id in {CHEBI_O2, CHEBI_H2O2}
            for participant in reaction.all_participants()
        )
        if has_oxygen_cochemistry:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen/peroxide chemistry indicates a different oxidoreductase branch",
            )

        produces_ammonia = any(
            participant.chebi_id in self.AMMONIA_CHEBIS
            for participant in reaction.right_participants
        )
        if produces_ammonia:
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia-releasing redox chemistry belongs to CH-NH2 donor classes",
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
        if not left_non_spectators or not right_non_spectators:
            return ClassificationResult(
                is_member=False,
                explanation="No non-cofactor substrate/product pair to evaluate",
            )

        has_oxo_donor = any(
            self._is_oxo_donor(participant) for participant in left_non_spectators
        )
        if not has_oxo_donor:
            return ClassificationResult(
                is_member=False,
                explanation="No aldehyde or oxo donor detected on the substrate side",
            )

        has_carboxylate_product = any(
            self._is_carboxylate_product(participant)
            for participant in right_non_spectators
        )
        has_thioester_product = any(
            participant.is_thioester() for participant in right_non_spectators
        )
        has_acyl_phosphate_product = any(
            self._is_acyl_phosphate_product(participant)
            for participant in right_non_spectators
        )
        has_co2_product = any(
            participant.chebi_id == CHEBI_CO2 for participant in right_non_spectators
        )

        if has_carboxylate_product or has_thioester_product or has_acyl_phosphate_product:
            return ClassificationResult(
                is_member=True,
                explanation="Aldehyde/oxo oxidoreductase: NAD(P)+-dependent oxidation to carboxylate, acyl phosphate, or thioester",
            )

        if any(self._is_formate_like(participant) for participant in left_non_spectators) and has_co2_product:
            return ClassificationResult(
                is_member=True,
                explanation="Aldehyde/oxo oxidoreductase: formate-like donor oxidized to CO2 with NAD(P)+ acceptor",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No oxidized product consistent with aldehyde/oxo donor chemistry",
        )

    @classmethod
    def _is_oxo_donor(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        if is_aldehyde(participant.smiles, participant.chebi_id) or cls._is_formate_like(participant):
            return True
        if not is_ketone(participant.smiles, participant.chebi_id):
            return False
        if participant.has_moiety(Moiety.HYDROXYL):
            return False
        return bool(
            is_carboxylic_acid(participant.smiles, participant.chebi_id)
            or participant.is_phosphorylated()
            or participant.is_thioester()
        )

    @staticmethod
    def _is_carboxylate_product(participant: Participant) -> bool:
        return bool(
            participant.smiles
            and is_carboxylic_acid(participant.smiles, participant.chebi_id)
        )

    @classmethod
    def _is_formate_like(cls, participant: Participant) -> bool:
        if not participant.smiles or cls.FORMATE_PATTERN is None:
            return False
        mol = participant.get_mol()
        return mol is not None and mol.HasSubstructMatch(cls.FORMATE_PATTERN)

    @classmethod
    def _is_acyl_phosphate_product(cls, participant: Participant) -> bool:
        if not participant.smiles or cls.ACYL_PHOSPHATE_PATTERN is None:
            return False
        mol = participant.get_mol()
        return mol is not None and mol.HasSubstructMatch(cls.ACYL_PHOSPHATE_PATTERN)
