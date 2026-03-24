"""CH-OH Oxidoreductase - EC 1.1 intermediate class.

EC 1.1: Oxidoreductases acting on the CH-OH group of donors.
These enzymes catalyze oxidation of alcohols to aldehydes/ketones.

This is the parent class for:
- Dehydrogenase (EC 1.1.1): NAD/NADP as acceptor
- AlcoholOxidase (EC 1.1.3): O2 as acceptor

Pattern: R-CH(OH)-R' + acceptor → R-C(=O)-R' + reduced acceptor

Classification strategy:
1. ChEBI IDs for exact compound matching (cofactors, ammonia)
2. SMARTS-based structural detection (aldehydes, enoyl, disulfides)
3. Mechanistic exclusions using cofactors and structural moieties
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_COA,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_NAD_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.moiety import (
    Moiety,
    has_disulfide,
    is_aldehyde,
    is_carboxylic_acid,
    is_enoyl,
)


class OxidoreductaseActingOnCHOHGroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on CH-OH group of donors

    These enzymes catalyze:
    - Oxidation: R-CHOH-R' → R-CO-R' (alcohol to carbonyl)
    - Reduction: R-CO-R' → R-CHOH-R' (carbonyl to alcohol)

    Key distinction from other oxidoreductase subclasses:
    - EC 1.1: CH-OH group (alcohol/hydroxyl)
    - EC 1.2: Aldehyde/oxo group (already oxidized carbon)
    - EC 1.3: CH-CH group (C=C bond formation/reduction)
    - EC 1.4: CH-NH2 group (amino group, releases ammonia)
    """

    GO_ID = "GO:0016614"  # oxidoreductase activity, acting on CH-OH group of donors
    EC_NUMBER_PREFIX = "1.1.-.-"

    PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838"}
    ACP_CHEBI_IDS = {"CHEBI:78784", "CHEBI:78785"}  # enoyl-ACP, acyl-ACP

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction acts on CH-OH group.

        Strategy:
        1. Must be an oxidoreductase (parent)
        2. Exclude CH-CH patterns (enoyl, unsaturated → EC 1.3)
        3. Exclude CH-NH2 patterns (ammonia release → EC 1.4)
        4. Exclude disulfide/thiol and ACP-linked chemistry
        5. Exclude oxygenase-like coupled redox (O2 + NAD(P)H + H2O/H2O2)
        6. Require NAD(P)+/NAD(P)H system
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Must use NAD+/NADP+ system for EC 1.1.1
        nad_system = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}
        has_nad_system = any(
            p.chebi_id in nad_system
            for p in reaction.left_participants + reaction.right_participants
        )

        if not has_nad_system:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD+/NADP+ system - EC 1.1.1 requires NAD/NADP"
            )

        # Exclude CH-CH patterns (these are EC 1.3, not EC 1.1)
        # Primary check: SMARTS-based enoyl detection
        has_enoyl_structure = any(
            p.smiles and is_enoyl(p.smiles, p.chebi_id)
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_enoyl_structure:
            return ClassificationResult(
                is_member=False,
                explanation="Enoyl group detected (SMARTS) - EC 1.3, not EC 1.1"
            )

        # Exclude CH-NH2 patterns (ammonia release → EC 1.4)
        ammonia_ids = {CHEBI_NH3, CHEBI_NH4}
        has_ammonia_by_id = any(
            p.chebi_id in ammonia_ids
            for p in reaction.right_participants
        )

        if has_ammonia_by_id:
            return ClassificationResult(
                is_member=False,
                explanation="Ammonia product (ChEBI) - EC 1.4 (CH-NH2), not EC 1.1"
            )

        # Exclude disulfide/thiol reactions
        has_disulfide_structure = any(
            p.smiles and has_disulfide(p.smiles, p.chebi_id)
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_disulfide_structure:
            return ClassificationResult(
                is_member=False,
                explanation="Disulfide bond detected (SMARTS) - not CH-OH"
            )

        # Exclude ACP-linked fatty-acid chemistry using ChEBI IDs or ACP SMARTS.
        has_acp = any(
            p.chebi_id in self.ACP_CHEBI_IDS or p.has_moiety(Moiety.ACP_LINKAGE)
            for p in reaction.left_participants + reaction.right_participants
        )
        if has_acp:
            return ClassificationResult(
                is_member=False,
                explanation="ACP-linked acyl chemistry - not CH-OH donor class"
            )

        # Exclude oxygenase-like coupled redox chemistry.
        # CH-OH oxidoreductases in this branch should not require O2 co-substrate
        # together with NAD(P)H and peroxide/water oxygen reduction products.
        has_o2_left = any(p.chebi_id == CHEBI_O2 for p in reaction.left_participants)
        has_reduced_nad_left = any(
            p.chebi_id in {CHEBI_NADH, CHEBI_NADPH} for p in reaction.left_participants
        )
        has_oxygen_reduction_products = any(
            p.chebi_id in {CHEBI_H2O, CHEBI_H2O2} for p in reaction.right_participants
        )
        if has_o2_left and has_reduced_nad_left and has_oxygen_reduction_products:
            return ClassificationResult(
                is_member=False,
                explanation="O2/NAD(P)H coupled oxygenase signature - not CH-OH donor redox"
            )

        # Exclude aldehyde/oxo oxidoreductase signatures (EC 1.2, not EC 1.1).
        # EC 1.2 is characterized by aldehyde <-> carboxylate interconversion.
        # In contrast, EC 1.1 alcohol dehydrogenases often produce aldehydes.
        has_aldehyde_left_structure = any(
            p.smiles and is_aldehyde(p.smiles, p.chebi_id)
            for p in reaction.left_participants
        )
        has_aldehyde_right_structure = any(
            p.smiles and is_aldehyde(p.smiles, p.chebi_id)
            for p in reaction.right_participants
        )
        has_carboxyl_left_structure = any(
            p.smiles and is_carboxylic_acid(p.smiles, p.chebi_id)
            for p in reaction.left_participants
        )
        has_carboxyl_right_structure = any(
            p.smiles and is_carboxylic_acid(p.smiles, p.chebi_id)
            for p in reaction.right_participants
        )

        if has_aldehyde_left_structure and has_carboxyl_right_structure:
            return ClassificationResult(
                is_member=False,
                explanation="Aldehyde -> carboxylate pattern (SMARTS) - EC 1.2, not EC 1.1"
            )
        if has_carboxyl_left_structure and has_aldehyde_right_structure:
            return ClassificationResult(
                is_member=False,
                explanation="Carboxylate -> aldehyde reverse EC 1.2 pattern (SMARTS)"
            )

        # Exclude acyl-phosphate aldehyde dehydrogenase chemistry (EC 1.2).
        # Mechanistic signature: aldehyde + Pi + NAD(P)+ -> acyl phosphate + NAD(P)H.
        has_phosphate_left = any(
            p.chebi_id in self.PHOSPHATE_IDS for p in reaction.left_participants
        )
        has_phosphate_right = any(
            p.chebi_id in self.PHOSPHATE_IDS for p in reaction.right_participants
        )
        if (has_phosphate_left and has_aldehyde_left_structure) or (
            has_phosphate_right and has_aldehyde_right_structure
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Aldehyde + phosphate signature (acyl-phosphate chemistry) - EC 1.2"
            )

        # Exclude CoA-coupled oxidative decarboxylation chemistry.
        has_coa = any(
            p.chebi_id == CHEBI_COA for p in reaction.left_participants + reaction.right_participants
        )
        has_co2 = any(
            p.chebi_id == CHEBI_CO2 for p in reaction.left_participants + reaction.right_participants
        )
        if has_coa and has_co2:
            return ClassificationResult(
                is_member=False,
                explanation="CoA-coupled decarboxylative redox chemistry - not CH-OH donor class"
            )

        # Exclude transhydrogenase (NAD+ + NADPH → NADH + NADP+)
        # This is just electron transfer between cofactors, not substrate oxidation
        nad_ids = {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}
        all_participants = reaction.left_participants + reaction.right_participants
        non_cofactor_participants = [
            p for p in all_participants
            if p.chebi_id not in nad_ids
        ]

        # If all participants are NAD/NADP cofactors, it's transhydrogenase
        if len(non_cofactor_participants) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P)+ transhydrogenase - electron transfer between cofactors"
            )

        # Exclude non-specific placeholders that are not chemically grounded.
        has_unresolved_participant = any(
            p.chebi_id is None and not p.smiles
            for p in non_cofactor_participants
        )
        if has_unresolved_participant:
            return ClassificationResult(
                is_member=False,
                explanation="Unresolved non-cofactor participants - cannot support CH-OH mechanistic call"
            )

        return ClassificationResult(
            is_member=True,
            explanation="CH-OH oxidoreductase: NAD(P)-dependent alcohol/carbonyl interconversion"
        )
