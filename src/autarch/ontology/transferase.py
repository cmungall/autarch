"""Simplified Transferase reaction classification.

Transferases (EC 2) catalyze the transfer of a functional group from one
molecule (donor) to another (acceptor). Key group types:
- Phosphoryl groups (kinases): ATP → ADP
- Methyl groups: SAM → SAH
- Glycosyl groups: NDP-sugar → NDP
- Acyl groups: acyl-CoA → CoA
- Amino groups: amino acid ↔ keto acid

History
-------
## 2025-12-27 (v4 - Simplified)
Major simplification to reduce complexity (CC 139 → target <30).
Uses declarative exclusion patterns and cleaner detection logic.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.molecules import (
    CHEBI_ATP, CHEBI_ADP, CHEBI_AMP, CHEBI_SAM, CHEBI_SAH, CHEBI_ACETYL_COA, CHEBI_COA,
    CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
    CHEBI_H2O, CHEBI_H_PLUS, CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, CHEBI_NH3, CHEBI_NH4, CHEBI_O2,
)

# Key ChEBI IDs for detection
CHEBI_BICARBONATE = "CHEBI:17544"  # hydrogencarbonate (HCO3-)
CHEBI_PAPS = "CHEBI:58339"  # 3'-phosphoadenylyl sulfate
CHEBI_PAP = "CHEBI:58343"   # adenosine 3',5'-bisphosphate

# Phosphate variants (ETL sometimes maps to polyphosphate CHEBI:16838)
PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838"}

# Amino/keto acid pairs for aminotransferase detection
KETO_ACIDS = {"CHEBI:16810", "CHEBI:16452", "CHEBI:15361"}  # 2-OG, OAA, pyruvate
AMINO_ACIDS = {"CHEBI:29985", "CHEBI:29991", "CHEBI:57972"}  # Glu, Asp, Ala

# NDP-sugar donors (glycosyltransferases)
NDP_SUGAR_DONORS = {
    "CHEBI:18066", "CHEBI:58885", "CHEBI:66914",  # UDP-glucose forms
    "CHEBI:18307", "CHEBI:67119",  # UDP-galactose
    "CHEBI:57705",  # UDP-glucuronate
    "CHEBI:17659", "CHEBI:57513",  # UDP-GlcNAc
    "CHEBI:17629", "CHEBI:57527",  # GDP-mannose
    "CHEBI:15819", "CHEBI:57273",  # GDP-fucose
    "CHEBI:16556",  # CMP-sialic acid
}

# NDP carriers (products of glycosyl transfer)
NDP_CARRIERS = {"CHEBI:58223", "CHEBI:17659", "CHEBI:18066"}  # UDP, GDP, etc


class Transferase(ReactionClass):
    """transferase"""

    GO_ID: ClassVar[Optional[str]] = "GO:0016740"
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "2.-.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for transferase activity using declarative patterns."""

        # === TRANSPORT EXCLUSION ===
        # Use the new location-aware transport detection
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not group transfer"
            )

        # === COLLECT PARTICIPANTS ===
        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        # === EXCLUSIONS (check these first) ===

        # 1. Carboxylase/ligase: HCO3 + ATP → product + ADP + Pi (carbon fixation)
        if CHEBI_BICARBONATE in left_ids and CHEBI_ATP in left_ids:
            if CHEBI_ADP in right_ids and (PHOSPHATE_IDS & right_ids):
                return ClassificationResult(
                    is_member=False,
                    explanation="Carboxylase/ligase - carbon fixation not group transfer"
                )

        # 2. Ligase: multiple substrates + ATP → product + AMP + PPi
        if CHEBI_ATP in left_ids and CHEBI_AMP in right_ids and CHEBI_DIPHOSPHATE in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Ligase - bond formation not group transfer"
            )

        # 2b. Ligase: multiple substrates + ATP → product + ADP + Pi (peptide bond formation)
        # These are ligases forming C-N or C-O bonds, not phosphoryl transfer
        if CHEBI_ATP in left_ids and CHEBI_ADP in right_ids and (PHOSPHATE_IDS & right_ids):
            # Count non-cofactor substrates (ligases join 2+ substrates)
            non_cofactor_left = left_ids - {CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS}
            if len(non_cofactor_left) >= 2:
                return ClassificationResult(
                    is_member=False,
                    explanation="Ligase - bond formation with multiple substrates, not phosphoryl transfer"
                )

        # 3. Oxidoreductase with NAD(P)+: excludes dehydrogenases misclassified
        has_nad = bool({CHEBI_NAD_PLUS, CHEBI_NADP_PLUS} & left_ids)
        has_nadh = bool({CHEBI_NADH, CHEBI_NADPH} & right_ids)
        if has_nad and has_nadh:
            return ClassificationResult(
                is_member=False,
                explanation="Oxidoreductase - NAD(P)+ dependent, not transferase"
            )

        # 4. Oxygenase: uses O2
        if CHEBI_O2 in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Oxygenase - uses O2, not transferase"
            )

        # 5. Deaminase: produces NH4+ (not aminotransferase)
        # Aminotransferases don't produce free ammonia - they transfer the amino group
        if {CHEBI_NH3, CHEBI_NH4} & right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Deaminase/lyase - produces ammonia, not amino group transfer"
            )

        # === POSITIVE DETECTION ===

        # 1. PHOSPHORYL TRANSFER (kinase): ATP → ADP (phosphate to substrate)
        if CHEBI_ATP in left_ids and CHEBI_ADP in right_ids:
            # Exclude hydrolysis: ATP + H2O → ADP + Pi (no other substrate)
            non_cofactor = left_ids - {CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS}
            if len(non_cofactor) >= 1:
                return ClassificationResult(
                    is_member=True,
                    explanation="Transferase: phosphoryl transfer (kinase)"
                )

        # 2. METHYL TRANSFER: SAM → SAH
        if CHEBI_SAM in left_ids and CHEBI_SAH in right_ids:
            return ClassificationResult(
                is_member=True,
                explanation="Transferase: methyl transfer (methyltransferase)"
            )

        # 3. SULFATE TRANSFER: PAPS + acceptor → PAP + sulfated-acceptor
        if CHEBI_PAPS in left_ids and CHEBI_PAP in right_ids:
            # Exclude hydrolysis: PAPS + H2O → sulfate + PAP (no acceptor)
            non_cofactor = left_ids - {CHEBI_PAPS, CHEBI_H2O, CHEBI_H_PLUS}
            if len(non_cofactor) >= 1:
                return ClassificationResult(
                    is_member=True,
                    explanation="Transferase: sulfate transfer (sulfotransferase)"
                )

        # 4. GLYCOSYL TRANSFER: NDP-sugar → NDP
        has_ndp_sugar = bool(NDP_SUGAR_DONORS & left_ids)
        # Check for NDP in products (or nucleotide pattern by SMILES)
        has_ndp_product = any(
            p.chebi_id in NDP_CARRIERS or
            (p.smiles and "O=C1C=CN" in p.smiles)  # Uracil pattern (UDP)
            for p in reaction.right_participants
        )
        if has_ndp_sugar and has_ndp_product:
            # Exclude epimerases: single NDP-sugar ↔ single NDP-sugar
            non_spectator_left = left_ids - {CHEBI_H_PLUS, CHEBI_H2O}
            non_spectator_right = right_ids - {CHEBI_H_PLUS, CHEBI_H2O}
            if len(non_spectator_left) == 1 and len(non_spectator_right) == 1:
                # Check if product is also NDP-sugar (epimerase)
                if NDP_SUGAR_DONORS & right_ids:
                    return ClassificationResult(
                        is_member=False,
                        explanation="Epimerase - NDP-sugar interconversion, not glycosyl transfer"
                    )
            return ClassificationResult(
                is_member=True,
                explanation="Transferase: glycosyl transfer (glycosyltransferase)"
            )

        # 5. ACYL TRANSFER: acyl-CoA + acceptor → CoA + acyl-acceptor
        has_coa_product = CHEBI_COA in right_ids
        has_thioester = any(p.is_thioester() for p in reaction.left_participants)
        has_acetyl_coa = CHEBI_ACETYL_COA in left_ids

        if has_coa_product or has_thioester or has_acetyl_coa:
            # Exclude single-substrate reactions (isomerases, lyases)
            non_cofactor_left_participants = [
                p for p in reaction.left_participants
                if p.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS, CHEBI_COA}
            ]
            non_cofactor_right_participants = [
                p for p in reaction.right_participants
                if p.chebi_id not in {CHEBI_H2O, CHEBI_H_PLUS, CHEBI_COA}
            ]

            # Single thioester → products is lyase/isomerase, not transfer
            if len(non_cofactor_left_participants) == 1 and has_thioester:
                # Check if it's just dehydration or isomerization
                if len(non_cofactor_right_participants) <= 2:
                    return ClassificationResult(
                        is_member=False,
                        explanation="CoA-dependent lyase/isomerase - not acyl transfer"
                    )

            # Need at least 2 substrates for true acyl transfer
            if len(non_cofactor_left_participants) >= 2:
                return ClassificationResult(
                    is_member=True,
                    explanation="Transferase: acyl transfer (acyltransferase)"
                )

        # 6. AMINOTRANSFERASE: amino acid + keto acid ↔ keto acid + amino acid
        # Must have reciprocal exchange (no NH4+ produced - checked above)
        left_has_keto = bool(KETO_ACIDS & left_ids)
        left_has_amino = bool(AMINO_ACIDS & left_ids)
        right_has_keto = bool(KETO_ACIDS & right_ids)
        right_has_amino = bool(AMINO_ACIDS & right_ids)

        # True aminotransferase: one amino in, different amino out (via keto intermediates)
        if (left_has_keto and right_has_amino) or (left_has_amino and right_has_keto):
            # Already excluded NH4+ production above
            return ClassificationResult(
                is_member=True,
                explanation="Transferase: amino group transfer (aminotransferase)"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No group transfer pattern detected"
        )
