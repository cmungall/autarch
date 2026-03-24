"""oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH group acts as hydrogen or electron donor and reduces oxygen.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import (
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
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_FERREDOXIN_OXIDIZED = "CHEBI:33737"
CHEBI_FERREDOXIN_REDUCED = "CHEBI:33738"
CHEBI_IRON_II = "CHEBI:29033"
CHEBI_IRON_III = "CHEBI:29034"


class OxidoreductaseActingOnTheCHCHGroupOfDonorsOxygenAsAcceptor(Oxidoreductase):
    """oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which a CH-CH group
    acts as hydrogen or electron donor and reduces oxygen.
    """

    GO_ID = "GO:0016634"
    EC_NUMBER_PREFIX = "1.3.3.-"
    EXCLUDED_CHEBIS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_FMN,
        CHEBI_FMNH2,
    }
    LEFT_SPECTATOR_CHEBIS = {CHEBI_O2, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H2O2, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for oxygen-acceptor CH-CH oxidation."""
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
        """Evaluate one reaction orientation."""
        if any(
            participant.chebi_id in self.EXCLUDED_CHEBIS
            for participant in left_participants + right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="External nicotinamide or flavin cofactors indicate a different CH-CH oxidoreductase branch",
            )

        if not any(participant.chebi_id == CHEBI_O2 for participant in left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No dioxygen acceptor in this reaction orientation",
            )

        if self._has_external_donor_pair(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="External iron-sulfur or cytochrome-like donor pairs indicate a different oxygenase branch",
            )

        if not any(
            participant.chebi_id in {CHEBI_H2O2, CHEBI_H2O}
            for participant in right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Oxygen acceptor branch requires peroxide or water products",
            )

        left_organics = self._organic_participants(
            left_participants,
            self.LEFT_SPECTATOR_CHEBIS,
        )
        right_organics = self._organic_participants(
            right_participants,
            self.RIGHT_SPECTATOR_CHEBIS,
        )
        if len(left_organics) != 1 or not right_organics:
            return ClassificationResult(
                is_member=False,
                explanation="Requires one substantive CH-CH donor substrate and structured products",
            )

        if self._total_oxygen(right_organics) > self._total_oxygen(left_organics):
            return ClassificationResult(
                is_member=False,
                explanation="Organic products gain oxygen overall, which is more consistent with oxygenase chemistry",
            )

        if self._count_unsaturated_cc_bonds(right_organics) <= self._count_unsaturated_cc_bonds(
            left_organics
        ):
            return ClassificationResult(
                is_member=False,
                explanation="No net increase in carbon-carbon unsaturation detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="CH-CH oxidoreductase with oxygen as acceptor",
        )

    @staticmethod
    def _has_external_donor_pair(
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        reduced_chebis = {CHEBI_FERREDOXIN_REDUCED, CHEBI_IRON_II}
        oxidized_chebis = {CHEBI_FERREDOXIN_OXIDIZED, CHEBI_IRON_III}
        has_reduced = any(participant.chebi_id in reduced_chebis for participant in left_participants)
        has_oxidized = any(
            participant.chebi_id in oxidized_chebis for participant in right_participants
        )
        return has_reduced and has_oxidized

    @staticmethod
    def _organic_participants(
        participants: list[Participant],
        excluded_chebis: set[str],
    ) -> list[Participant]:
        organics: list[Participant] = []
        for participant in participants:
            if participant.chebi_id in excluded_chebis:
                continue
            mol = participant.get_mol()
            if mol is None:
                continue
            if any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
                organics.append(participant)
        return organics

    @staticmethod
    def _total_oxygen(participants: list[Participant]) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            total += oxygen_atoms * max(participant.count, 1)
        return total

    @staticmethod
    def _count_unsaturated_cc_bonds(participants: list[Participant]) -> int:
        total = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            for bond in mol.GetBonds():
                begin = bond.GetBeginAtom().GetAtomicNum()
                end = bond.GetEndAtom().GetAtomicNum()
                if begin != 6 or end != 6:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetIsAromatic():
                    total += max(participant.count, 1)
        return total
