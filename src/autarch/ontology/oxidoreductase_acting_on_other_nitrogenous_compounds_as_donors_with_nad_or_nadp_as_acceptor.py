"""oxidoreductase acting on other nitrogenous compounds as donors with NAD or NADP as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a nitrogenous
 group, excluding NH and NH2 groups, acts as a hydrogen or electron donor and
 reduces NAD or NADP.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh2_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHNH2GroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_nh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor,
)


class OxidoreductaseActingOnOtherNitrogenousCompoundsAsDonorsWithNADOrNADPAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase acting on other nitrogenous compounds as donors with NAD or NADP as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a
    nitrogenous group, excluding NH and NH2 groups, acts as a hydrogen or
    electron donor and reduces NAD or NADP.
    """

    GO_ID = "GO:0046857"
    EC_NUMBER_PREFIX = "1.7.1.-"
    OXIDIZED_NICOTINAMIDES = {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
    REDUCED_NICOTINAMIDES = {CHEBI_NADH, CHEBI_NADPH}
    SPECTATOR_CHEBIS = (
    OXIDIZED_NICOTINAMIDES
        | REDUCED_NICOTINAMIDES
        | {CHEBI_H2O, CHEBI_H_PLUS}
    )
    INORGANIC_NITROGEN_OXIDES = {
        "CHEBI:16301",  # nitrite
        "CHEBI:17632",  # nitrate
        "CHEBI:15429",  # hydroxylamine
        "CHEBI:16480",  # nitric oxide
        "CHEBI:17045",  # nitrous oxide
    }
    N_O_PATTERN = Chem.MolFromSmarts("[N;!R]-[O]")
    N_OXIDE_PATTERN = Chem.MolFromSmarts("[n+][O-]")
    AZO_PATTERN = Chem.MolFromSmarts("[N]=[N]")
    AZOXY_PATTERN = Chem.MolFromSmarts("[N]-[N]=[O]")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for NAD(P)-linked oxidation of non-CH-N nitrogenous donors."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not a nitrogenous-compound oxidoreductase",
            )

        if self._has_external_oxidant(reaction):
            return ClassificationResult(
                is_member=False,
                explanation="External oxygen/peroxide oxidant indicates a different oxidoreductase branch",
            )

        if not self._has_nicotinamide_pair(reaction.left_participants, reaction.right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Requires NAD(+) or NADP(+) reduced across the reaction",
            )

        ch_nh = OxidoreductaseActingOnTheCHNHGroupOfDonorsNADOrNADPAsAcceptor()
        if ch_nh.check_membership_impl(reaction).is_member:
            return ClassificationResult(
                is_member=False,
                explanation="Classified more specifically as a CH-NH oxidoreductase",
            )

        ch_nh2 = OxidoreductaseActingOnTheCHNH2GroupOfDonorsNADOrNADPAsAcceptor()
        if ch_nh2.check_membership_impl(reaction).is_member:
            return ClassificationResult(
                is_member=False,
                explanation="Classified more specifically as a CH-NH2 oxidoreductase",
            )

        reactive = [
            participant
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id not in self.SPECTATOR_CHEBIS
        ]
        nitrogenous = [participant for participant in reactive if self._nitrogen_count(participant) > 0]
        if not nitrogenous:
            return ClassificationResult(
                is_member=False,
                explanation="No substantive nitrogen-containing participants remain after removing cofactors",
            )

        if self._has_other_nitrogenous_redox_signature(nitrogenous):
            return ClassificationResult(
                is_member=True,
                explanation="NAD(P)-linked oxidoreduction of a non-CH-N/NH2 nitrogenous donor",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No qualifying non-CH-N nitrogenous redox signature detected",
        )

    @classmethod
    def _has_nicotinamide_pair(
        cls,
        left: list[Participant],
        right: list[Participant],
    ) -> bool:
        left_oxidized = cls._count_chebis(left, cls.OXIDIZED_NICOTINAMIDES)
        right_reduced = cls._count_chebis(right, cls.REDUCED_NICOTINAMIDES)
        left_reduced = cls._count_chebis(left, cls.REDUCED_NICOTINAMIDES)
        right_oxidized = cls._count_chebis(right, cls.OXIDIZED_NICOTINAMIDES)
        return left_oxidized in {1, 2, 3} and left_oxidized == right_reduced and left_reduced == 0 and right_oxidized == 0

    @staticmethod
    def _count_chebis(participants: list[Participant], chebis: set[str]) -> int:
        return sum(max(participant.count, 1) for participant in participants if participant.chebi_id in chebis)

    @staticmethod
    def _has_external_oxidant(reaction: Reaction) -> bool:
        chebis = {
            participant.chebi_id
            for participant in reaction.left_participants + reaction.right_participants
            if participant.chebi_id
        }
        return CHEBI_O2 in chebis or CHEBI_H2O2 in chebis

    @classmethod
    def _has_other_nitrogenous_redox_signature(cls, participants: list[Participant]) -> bool:
        for participant in participants:
            if participant.chebi_id in cls.INORGANIC_NITROGEN_OXIDES:
                return True
            if cls._has_pattern(participant, cls.N_O_PATTERN):
                return True
            if cls._has_pattern(participant, cls.N_OXIDE_PATTERN):
                return True
            if cls._has_pattern(participant, cls.AZO_PATTERN):
                return True
            if cls._has_pattern(participant, cls.AZOXY_PATTERN):
                return True
        return False

    @staticmethod
    def _nitrogen_count(participant: Participant) -> int:
        mol = participant.get_mol()
        if mol is None:
            return 0
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    @staticmethod
    def _has_pattern(participant: Participant, pattern: Chem.Mol | None) -> bool:
        if pattern is None:
            return False
        mol = participant.get_mol()
        return mol is not None and mol.HasSubstructMatch(pattern)
