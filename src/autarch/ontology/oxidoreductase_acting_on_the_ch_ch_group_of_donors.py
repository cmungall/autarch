"""CH-CH Oxidoreductase - EC 1.3 intermediate class.

EC 1.3: Oxidoreductases acting on the CH-CH group of donors.
These enzymes catalyze saturation/desaturation of C-C bonds.

This is the parent class for:
- EnoylReductase (EC 1.3.1): NAD/NADP as acceptor
- Desaturase (EC 1.3.8): FAD as acceptor

Pattern: R-CH=CH-R' + acceptor-H2 → R-CH2-CH2-R' + acceptor
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_O2,
)
from autarch.moiety import Moiety, is_aldehyde, is_enoyl


class OxidoreductaseActingOnTheCHCHGroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on the CH-CH group of donors

    These enzymes catalyze:
    - Desaturation: R-CH2-CH2-R' → R-CH=CH-R' (single to double bond)
    - Saturation: R-CH=CH-R' → R-CH2-CH2-R' (double to single bond)

    Key distinction from other oxidoreductase subclasses:
    - EC 1.1: CH-OH group (alcohol oxidation/reduction)
    - EC 1.2: Aldehyde/oxo group
    - EC 1.3: CH-CH group (C=C bond interconversion) ← THIS CLASS
    - EC 1.4: CH-NH2 group (amino group)

    Characteristic patterns:
    - Enoyl-CoA/acyl-CoA interconversion (fatty acid metabolism)
    - Dihydro-compounds to aromatic/unsaturated (dihydroorotate → orotate)
    - Saturated → unsaturated naming patterns
    """

    GO_ID = "GO:0016627"  # oxidoreductase activity, acting on the CH-CH group of donors
    EC_NUMBER_PREFIX = "1.3.-.-"
    KETONE_PATTERN = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    SPECTATOR_IDS = {
        CHEBI_NAD_PLUS,
        CHEBI_NADH,
        CHEBI_NADP_PLUS,
        CHEBI_NADPH,
        CHEBI_FAD,
        CHEBI_FADH2,
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction acts on CH-CH group.

        Strategy:
        1. Must be an oxidoreductase (parent)
        2. Require redox cofactor system (NAD(P)/FAD)
        3. Detect CH-CH chemistry structurally:
           - enoyl motif, or
           - change in C=C bond count across substrate/product sides
        4. Exclude pure CH-OH (alcohol <-> carbonyl) signatures without CH-CH evidence
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Must use NAD+/NADP+ or FAD system
        nad_system = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}
        fad_system = {CHEBI_FAD, CHEBI_FADH2}
        redox_cofactors = nad_system | fad_system

        has_redox_system = any(
            p.chebi_id in redox_cofactors
            for p in reaction.left_participants + reaction.right_participants
        )

        if not has_redox_system:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD/NADP/FAD system - EC 1.3 requires these cofactors"
            )

        # Exclude oxygenase-like coupled chemistry:
        # O2 + reduced redox cofactor -> H2O/H2O2 is typically oxygenase/monooxygenase.
        has_o2_left = any(p.chebi_id == CHEBI_O2 for p in reaction.left_participants)
        has_reduced_left = any(
            p.chebi_id in {CHEBI_NADH, CHEBI_NADPH, CHEBI_FADH2}
            for p in reaction.left_participants
        )
        has_oxidized_right = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS, CHEBI_FAD}
            for p in reaction.right_participants
        )
        has_oxygen_reduction_products = any(
            p.chebi_id in {CHEBI_H2O, CHEBI_H2O2}
            for p in reaction.right_participants
        )
        if has_o2_left and has_reduced_left and (has_oxidized_right or has_oxygen_reduction_products):
            return ClassificationResult(
                is_member=False,
                explanation="O2-coupled oxygenase signature - not CH-CH donor oxidoreductase"
            )

        left_non_spectators = [
            p for p in reaction.left_participants if p.chebi_id not in self.SPECTATOR_IDS
        ]
        right_non_spectators = [
            p for p in reaction.right_participants if p.chebi_id not in self.SPECTATOR_IDS
        ]
        all_non_spectators = left_non_spectators + right_non_spectators

        # Structural CH-CH evidence:
        # 1) alpha-beta unsaturated acyl motif (enoyl), or
        # 2) net change in C=C bond counts across sides.
        has_enoyl_structure = any(
            p.smiles and is_enoyl(p.smiles, p.chebi_id)
            for p in all_non_spectators
        )
        left_cc_double_bonds = self._count_cc_double_bonds(left_non_spectators)
        right_cc_double_bonds = self._count_cc_double_bonds(right_non_spectators)
        has_cc_interconversion = left_cc_double_bonds != right_cc_double_bonds

        if not (has_enoyl_structure or has_cc_interconversion):
            return ClassificationResult(
                is_member=False,
                explanation="No structural CH-CH interconversion signature detected"
            )

        # CH-OH exclusion: if the non-spectator conversion is purely alcohol<->carbonyl
        # and we lack independent CH-CH evidence, route to EC 1.1.
        has_alcohol_left = any(self._is_alcohol_like(p) for p in left_non_spectators)
        has_alcohol_right = any(self._is_alcohol_like(p) for p in right_non_spectators)
        has_carbonyl_left = any(self._is_carbonyl_like(p) for p in left_non_spectators)
        has_carbonyl_right = any(self._is_carbonyl_like(p) for p in right_non_spectators)
        has_choh_oxidoreduction = (has_alcohol_left and has_carbonyl_right) or (
            has_carbonyl_left and has_alcohol_right
        )
        if has_choh_oxidoreduction and not has_cc_interconversion and not has_enoyl_structure:
            return ClassificationResult(
                is_member=False,
                explanation="Alcohol/carbonyl redox without CH-CH bond change - EC 1.1, not EC 1.3"
            )

        # Determine acceptor type
        has_fad = any(
            p.chebi_id in fad_system
            for p in reaction.left_participants + reaction.right_participants
        )

        acceptor = "FAD" if has_fad else "NAD(P)+"

        return ClassificationResult(
            is_member=True,
            explanation=f"CH-CH oxidoreductase: {acceptor}-dependent saturation/desaturation"
        )

    @classmethod
    def _count_cc_double_bonds(cls, participants: list[Participant]) -> int:
        """Count explicit carbon-carbon double bonds across participants."""
        count = 0
        for participant in participants:
            mol = participant.get_mol()
            if mol is None:
                continue
            for bond in mol.GetBonds():
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                    count += 1
        return count

    @classmethod
    def _is_alcohol_like(cls, participant: Participant) -> bool:
        return participant.has_moiety(Moiety.HYDROXYL)

    @classmethod
    def _is_carbonyl_like(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        if is_aldehyde(participant.smiles, participant.chebi_id):
            return True
        mol = participant.get_mol()
        if mol is None:
            return False
        return cls.KETONE_PATTERN is not None and mol.HasSubstructMatch(cls.KETONE_PATTERN)
