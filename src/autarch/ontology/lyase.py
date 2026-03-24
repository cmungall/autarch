"""Simplified Lyase reaction classification.

Lyases cleave bonds without hydrolysis or oxidation, often forming double bonds.

History
-------

## 2025-12-23 (v3)

Fixed false positives by:
1. Refined stoichiometry rule: lyases typically 1→2 or 1→3 (one substrate splits)
   NOT 2→3 which is often transferase + leaving group
2. Added PAPS/PAP (sulfotransferase cofactors) to exclusion list
3. Added phosphate exclusion for phosphorolytic reactions

## 2025-12-22

Fixed false positives by:
1. Using correct ChEBI IDs for CoA cofactors (CHEBI:57287, CHEBI:57288)
2. Adding UDP/GDP nucleotide exclusion for glycosyltransferases
3. Being more strict about the "molecule count increased" fallback
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass

from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_SAH,
    CHEBI_SAM,
    CHEBI_UDP,
    CHEBI_GDP,
)

# Ferredoxin redox cofactors
CHEBI_FERREDOXIN_OX = "CHEBI:33737"  # oxidized [2Fe-2S]-ferredoxin
CHEBI_FERREDOXIN_RED = "CHEBI:33738"  # reduced [2Fe-2S]-ferredoxin

# Dihydrogen for hydrogenase detection
CHEBI_H2 = "CHEBI:18276"

# Malonyl-CoA for polyketide synthase detection
CHEBI_MALONYL_COA = "CHEBI:57384"

# Sulfotransferase cofactors
CHEBI_PAPS = "CHEBI:58339"  # 3'-phosphoadenylyl sulfate
CHEBI_PAP = "CHEBI:58343"   # adenosine 3',5'-bisphosphate

# Correct CoA ChEBI IDs used in RHEA
CHEBI_COA = "CHEBI:57287"  # coenzyme A(4-)
CHEBI_ACETYL_COA = "CHEBI:57288"  # acetyl-CoA(4-)


class Lyase(ReactionClass):
    """lyase"""
    
    GO_ID = "GO:0016829"  # lyase activity
    EC_NUMBER_PREFIX = "4.-.-.-"  # All lyases
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a lyase.
        
        Simple criteria: 
        - CO2 release (decarboxylation)
        - Water elimination (dehydration)
        - No water as reactant (distinguishes from hydrolase)
        """
        # Check for water as reactant (excludes hydrolases)
        has_water_reactant = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )
        
        if has_water_reactant:
            return ClassificationResult(
                is_member=False,
                explanation="Water as reactant - likely hydrolase not lyase"
            )
        
        # Check for CO2 as product (decarboxylation)
        has_co2_product = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.right_participants
        )
        
        if has_co2_product:
            # Check it's not oxidative decarboxylation (would be oxidoreductase)
            has_redox_cofactor = any(
                p.chebi_id in {
                    CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
                    CHEBI_FERREDOXIN_OX, CHEBI_FERREDOXIN_RED,  # ferredoxin
                    CHEBI_H2,  # dihydrogen (hydrogenase)
                }
                for p in reaction.left_participants + reaction.right_participants
            )

            # Also check for oxygen consumption (indicates oxidative process)
            has_oxygen = any(
                p.chebi_id == CHEBI_O2
                for p in reaction.left_participants
            )

            # Exclude polyketide synthases (malonyl-CoA + CO2 release)
            has_malonyl_coa = any(
                p.chebi_id == CHEBI_MALONYL_COA
                for p in reaction.left_participants
            )

            # Exclude CoA-dependent synthases (CoA in products with CO2)
            has_coa_product = any(
                p.chebi_id == CHEBI_COA
                for p in reaction.right_participants
            )

            # Exclude 2-oxoglutarate dehydrogenase complex (oxidative decarboxylation)
            CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
            has_2_oxoglutarate = any(
                p.chebi_id == CHEBI_2_OXOGLUTARATE
                for p in reaction.left_participants + reaction.right_participants
            )

            if not has_redox_cofactor and not has_oxygen and not has_malonyl_coa and not has_coa_product and not has_2_oxoglutarate:
                return ClassificationResult(
                    is_member=True,
                    explanation="Lyase: non-oxidative decarboxylation (CO2 release)"
                )
        
        # Check for water as product (dehydration)
        has_water_product = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )
        
        if has_water_product:
            # Check it's not an oxidase reaction (oxidases often produce water)
            has_redox_cofactor = any(
                p.chebi_id in {
                    CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
                    CHEBI_FERREDOXIN_OX, CHEBI_FERREDOXIN_RED, CHEBI_H2,
                }
                for p in reaction.left_participants + reaction.right_participants
            )
            
            has_oxygen = any(
                p.chebi_id == CHEBI_O2
                for p in reaction.left_participants
            )
            
            # Also check for transferase/ligase cofactors
            has_transferase_cofactor = any(
                p.chebi_id in {
                    CHEBI_SAM, CHEBI_SAH,  # methyltransferase cofactors
                    CHEBI_COA, CHEBI_ACETYL_COA,  # acyltransferase cofactors
                    CHEBI_UDP, CHEBI_GDP,  # glycosyltransferase cofactors
                    CHEBI_ATP, CHEBI_ADP,  # ligase cofactors (ATP-dependent synthesis)
                }
                for p in reaction.left_participants + reaction.right_participants
            )

            # Exclude 2-oxoglutarate dependent reactions (aminotransferases)
            CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
            has_2_oxoglutarate = any(
                p.chebi_id == CHEBI_2_OXOGLUTARATE
                for p in reaction.left_participants + reaction.right_participants
            )

            if not has_redox_cofactor and not has_oxygen and not has_transferase_cofactor and not has_2_oxoglutarate:
                return ClassificationResult(
                    is_member=True,
                    explanation="Lyase: dehydration (water elimination)"
                )
        
        # Check for ammonia release (common in C-N lyases)
        has_ammonia_product = any(
            p.chebi_id in {CHEBI_NH3, CHEBI_NH4}  # NH3, NH4+
            for p in reaction.right_participants
        )

        if has_ammonia_product and not has_water_reactant:
            # Check it's not an oxidative deamination (would be oxidoreductase)
            has_redox_cofactor = any(
                p.chebi_id in {
                    CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
                    CHEBI_FERREDOXIN_OX, CHEBI_FERREDOXIN_RED, CHEBI_H2,
                }
                for p in reaction.left_participants + reaction.right_participants
            )

            # Exclude 2-oxoglutarate dependent (aminotransferases)
            CHEBI_2_OXOGLUTARATE = "CHEBI:16810"
            has_2_oxoglutarate = any(
                p.chebi_id == CHEBI_2_OXOGLUTARATE
                for p in reaction.left_participants + reaction.right_participants
            )

            if not has_redox_cofactor and not has_2_oxoglutarate:
                return ClassificationResult(
                    is_member=True,
                    explanation="Lyase: C-N bond cleavage (ammonia release)"
                )

        # PRINCIPLED STOICHIOMETRY CHECK:
        # Lyases characteristically: one complex molecule ↔ multiple simpler molecules
        # - 1 → 2 or 1 → 3 (splitting)
        # - 2 → 1 or 3 → 1 (addition)
        # NOT 2 → 3 which is typically transferase + leaving group
        left_count = len(reaction.left_participants)
        right_count = len(reaction.right_participants)

        # Lyase-like stoichiometry: one side has 1 molecule, other has 2+
        is_lyase_stoichiometry = (
            (left_count == 1 and right_count >= 2) or  # splitting
            (left_count >= 2 and right_count == 1)     # addition
        )

        if is_lyase_stoichiometry:
            # SMILES patterns for nucleotide bases (many RHEA entries lack ChEBI IDs)
            UDP_SMILES_PATTERN = "O=C1C=CN"  # Uracil base pattern
            GDP_SMILES_PATTERN = "NC1=NC2=C(N=CN2"  # Guanine base pattern
            CMP_SMILES_PATTERN = "NC1=NC(=O)N"  # Cytosine base pattern

            def has_nucleotide_smiles(participant) -> bool:
                """Check if participant SMILES contains nucleotide pattern."""
                if not participant.smiles:
                    return False
                return (UDP_SMILES_PATTERN in participant.smiles or
                        GDP_SMILES_PATTERN in participant.smiles or
                        CMP_SMILES_PATTERN in participant.smiles)

            # Exclude reactions with redox cofactors (oxidoreductases)
            has_redox_cofactor = any(
                p.chebi_id in {
                    CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
                    CHEBI_FERREDOXIN_OX, CHEBI_FERREDOXIN_RED, CHEBI_H2,
                }
                for p in reaction.left_participants + reaction.right_participants
            )

            # Exclude reactions with transferase cofactors (check both ChEBI IDs and SMILES)
            has_transferase_cofactor = any(
                p.chebi_id in {
                    CHEBI_SAM, CHEBI_SAH,  # methyltransferase cofactors
                    CHEBI_COA, CHEBI_ACETYL_COA,  # acyltransferase cofactors
                    CHEBI_UDP, CHEBI_GDP,  # glycosyltransferase cofactors
                    CHEBI_PAPS, CHEBI_PAP,  # sulfotransferase cofactors
                    CHEBI_ATP,  # ligase/kinase cofactor
                    CHEBI_PHOSPHATE,  # phosphorolysis (not lyase)
                } or has_nucleotide_smiles(p)  # SMILES-based detection for nucleotides
                for p in reaction.left_participants + reaction.right_participants
            )

            # Exclude oxidative reactions
            has_oxygen = any(
                p.chebi_id == CHEBI_O2
                for p in reaction.left_participants
            )

            if not has_water_reactant and not has_redox_cofactor and not has_transferase_cofactor and not has_oxygen:
                return ClassificationResult(
                    is_member=True,
                    explanation="Lyase: non-hydrolytic bond cleavage/formation"
                )
        
        return ClassificationResult(
            is_member=False,
            explanation="No lyase pattern detected"
        )
