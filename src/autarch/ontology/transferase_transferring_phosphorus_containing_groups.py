"""Phosphorus-containing group transferase classification.

This broader class encompasses all EC 2.7 enzymes that transfer phosphorus-containing groups,
including kinases, nucleotidyltransferases, and other phosphotransferases.

History
-------

## 2025-12-22

Fixed false positives by:
1. Enabling ligase exclusion (ATP + substrates → product + AMP + PPi)
2. Enabling transport ATPase exclusion (ion movement across membranes)
3. Adding aminoacyl-tRNA synthetase detection
4. Calling parent Transferase class for proper validation
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase

from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CO2,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GMP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_PHOSPHATE,
    CHEBI_UTP,
    CHEBI_UDP,
)

# Bicarbonate for carboxylase exclusion
CHEBI_BICARBONATE = "CHEBI:17544"
PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838"}
DIPHOSPHATE_IDS = {CHEBI_DIPHOSPHATE, "CHEBI:18036"}


class TransferaseTransferringPhosphorusContainingGroups(Transferase):
    """transferase transferring phosphorus-containing groups
    
    This includes:
    - EC 2.7.1 - Phosphotransferases with alcohol group as acceptor (kinases)
    - EC 2.7.2 - Phosphotransferases with carboxy group as acceptor  
    - EC 2.7.3 - Phosphotransferases with nitrogenous group as acceptor
    - EC 2.7.4 - Phosphotransferases with phosphate group as acceptor
    - EC 2.7.7 - Nucleotidyltransferases (produce diphosphate)
    - EC 2.7.9 - Phosphotransferases with paired acceptors
    - EC 2.7.10-14 - Protein kinases
    """
    
    GO_ID = "GO:0016772"  # transferase activity, transferring phosphorus-containing groups
    EC_NUMBER_PREFIX = "2.7.-.-"  # All phosphotransferases
    
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction involves phosphorus-containing group transfer.
        
        This is a broad class that includes:
        1. Kinases (ATP/GTP → ADP/GDP)
        2. Nucleotidyltransferases (NTP → PPi)
        3. Other phosphotransferases
        """
        
        # Common phosphorus-containing donors
        PHOSPHO_DONORS = {
            CHEBI_ATP,  # ATP
            CHEBI_GTP,  # GTP
            CHEBI_CTP,  # CTP
            CHEBI_UTP,  # UTP
            CHEBI_ADP,  # ADP (for reverse reactions)
            CHEBI_GDP,   # GDP
            CHEBI_CDP,   # CDP
            CHEBI_UDP,   # UDP
            *PHOSPHATE_IDS,   # phosphate (for some transferases)
            *DIPHOSPHATE_IDS,   # diphosphate (PPi)
        }
        
        # Common phosphorus-containing products
        PHOSPHO_PRODUCTS = {
            CHEBI_ADP,  # ADP
            CHEBI_GDP,   # GDP
            CHEBI_CDP,   # CDP
            CHEBI_UDP,   # UDP
            CHEBI_ATP,   # ATP (for reverse reactions)
            CHEBI_GTP,   # GTP
            CHEBI_CTP,   # CTP
            CHEBI_UTP,   # UTP
            *DIPHOSPHATE_IDS,   # diphosphate (PPi)
            *PHOSPHATE_IDS,   # phosphate
            CHEBI_AMP,  # AMP
            CHEBI_GMP,   # GMP
        }
        
        # Check for phosphorus donors in reactants
        has_phospho_donor = any(
            p.chebi_id in PHOSPHO_DONORS for p in reaction.left_participants
        )
        
        # Check for phosphorus products
        has_phospho_product = any(
            p.chebi_id in PHOSPHO_PRODUCTS for p in reaction.right_participants
        )
        
        if has_phospho_donor and has_phospho_product:
            # Apply exclusions before classifying as phosphotransferase
            
            # Exclude hydrolysis reactions (phosphatases) - use water to cleave phosphate bonds
            has_water = any(
                p.chebi_id == CHEBI_H2O
                for p in reaction.left_participants
            )
            
            # If water is present and we're producing phosphate (not transferring), it's likely a phosphatase
            produces_free_phosphate = any(
                p.chebi_id in PHOSPHATE_IDS
                for p in reaction.right_participants
            )
            
            if has_water and produces_free_phosphate:
                # Check if this is simple phosphate hydrolysis rather than transfer
                non_water_reactants = [
                    p for p in reaction.left_participants
                    if p.chebi_id != CHEBI_H2O]
                
                if len(non_water_reactants) == 1:
                    return ClassificationResult(
                        is_member=False,
                        explanation="Phosphatase reaction - hydrolysis not transfer"
                    )
            
            # Define ATP-related variables early for use in exclusions
            has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
            has_amp = any(p.chebi_id == CHEBI_AMP for p in reaction.right_participants)
            has_ppi = any(
                p.chebi_id in DIPHOSPHATE_IDS for p in reaction.right_participants
            )
            has_adp = any(p.chebi_id == CHEBI_ADP for p in reaction.right_participants)

            # Exclude ligase reactions (bond formation using ATP)
            # Pattern: ATP + substrate(s) → product + AMP + PPi
            # ATP → AMP + PPi is the classic ligase signature (EC 6.x.x.x)
            if has_atp and has_amp and has_ppi:
                return ClassificationResult(
                    is_member=False,
                    explanation="Ligase reaction - ATP→AMP+PPi indicates bond formation not phosphoryl transfer"
                )

            # Exclude carboxylases (HCO3 + ATP → product + ADP + Pi)
            has_bicarbonate = any(
                p.chebi_id == CHEBI_BICARBONATE
                for p in reaction.left_participants
            )
            if has_bicarbonate and has_atp:
                return ClassificationResult(
                    is_member=False,
                    explanation="Carboxylase - ATP-dependent CO2 fixation, not phosphoryl transfer"
                )

            # Exclude PEPCK-like decarboxylases (oxaloacetate + GTP → PEP + GDP + CO2)
            # These use GTP for phosphoryl transfer to create PEP, but are lyases/decarboxylases
            has_gtp = any(p.chebi_id == CHEBI_GTP for p in reaction.left_participants)
            has_gdp = any(p.chebi_id == CHEBI_GDP for p in reaction.right_participants)
            has_co2 = any(
                p.chebi_id == CHEBI_CO2
                for p in reaction.right_participants
            )
            if has_gtp and has_gdp and has_co2:
                return ClassificationResult(
                    is_member=False,
                    explanation="PEPCK/decarboxylase - GTP-dependent decarboxylation, not simple phosphotransferase"
                )

            # CRITICAL: True phosphotransferases don't release free phosphate!
            # They transfer the phosphate TO another molecule.
            # If free Pi is released with ATP→ADP, it's:
            # - ATPase (hydrolysis): ATP + H2O → ADP + Pi
            # - Ligase: substrates + ATP → joined product + ADP + Pi
            # - Amidase/hydrolase with ATP activation
            if has_atp and has_adp and produces_free_phosphate:
                # Count non-ATP, non-water substrates
                non_cofactor = [
                    p for p in reaction.left_participants
                    if p.chebi_id not in {CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS}
                ]

                # If there are other substrates + ATP + H2O → products + ADP + Pi,
                # this is likely ATP-powered hydrolysis/ligase, not phosphotransferase
                if has_water:
                    return ClassificationResult(
                        is_member=False,
                        explanation="ATP hydrolysis with other products - ligase/amidase not phosphotransferase"
                    )

                # Even without water, if free Pi is released with substrates being transformed,
                # it's likely ligase (bond formation), not kinase (phosphoryl transfer)
                if len(non_cofactor) >= 1:
                    return ClassificationResult(
                        is_member=False,
                        explanation="Ligase - ATP powers bond formation, releasing free phosphate"
                    )

            # Use location-aware transport detection
            if reaction.is_transport_reaction():
                return ClassificationResult(
                    is_member=False,
                    explanation="Transport ATPase - moves molecules across membrane"
                )

            # Fallback: Check for shared molecules (transport pattern)
            left_chebi = {pr.chebi_id for pr in reaction.left_participants if pr.chebi_id}
            right_chebi = {pr.chebi_id for pr in reaction.right_participants if pr.chebi_id}
            spectator_ids = {CHEBI_ATP, CHEBI_ADP, CHEBI_H2O, CHEBI_H_PLUS, *PHOSPHATE_IDS}
            shared_molecules = (left_chebi & right_chebi) - spectator_ids

            if shared_molecules and has_atp and has_water:
                return ClassificationResult(
                    is_member=False,
                    explanation="Transport ATPase - moves molecules not transfers phosphorus"
                )

            # Exclude oxidoreductase reactions (use NADP+/NADPH, O2, etc.)
            has_redox_cofactor = any(
                p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH} for p in reaction.left_participants + reaction.right_participants
            )

            if has_redox_cofactor:
                return ClassificationResult(
                    is_member=False,
                    explanation="Oxidoreductase reaction - not phosphoryl transfer"
                )

            # Now classify as phosphotransferase
            explanation = "Phosphotransferase: "
            has_gtp = any(p.chebi_id == CHEBI_GTP for p in reaction.left_participants)
            has_gdp = any(p.chebi_id == CHEBI_GDP for p in reaction.right_participants)
            
            if has_ppi:
                explanation += "nucleotidyltransferase (produces diphosphate)"
            elif has_atp and has_adp:
                explanation += "ATP-dependent phosphoryl transfer"
            elif has_gtp and has_gdp:
                explanation += "GTP-dependent phosphoryl transfer"
            else:
                explanation += "phosphorus-containing group transfer"
            
            return ClassificationResult(
                is_member=True,
                explanation=explanation
            )

        return ClassificationResult(
            is_member=False,
            explanation="No phosphorus-containing group transfer detected"
        )
