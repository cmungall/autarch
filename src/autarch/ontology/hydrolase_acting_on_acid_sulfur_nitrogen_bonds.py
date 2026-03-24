"""Hydrolase acting on sulfur-nitrogen bonds (GO:0016826, EC 3.10.-.-) using pattern DSL.

This is a specific subclass of hydrolase that catalyzes the hydrolysis of
acid sulfur-nitrogen bonds.

History
-------

## 2025-12-22

Fixed false positives by:
1. Checking O-S sulfate ester patterns BEFORE inferring S-N from elements
2. The pattern "N-acetyl-D-galactosamine 4-sulfate" has S, N, produces sulfate,
   AND has "amine" in name - but is O-S hydrolysis, not S-N

## 2025-12-21

Refactored to use SMARTS-based sulfamate detection (Moiety.SULFAMATE) via
Participant.is_sulfamate() for structural S-N bond recognition. Name-based
patterns retained as fallback for molecules without SMILES data.
"""

from abc import abstractmethod

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import water, h_plus, p


class HydrolaseActingOnAcidSulfurNitrogenBonds(Hydrolase):
    """Legacy implementation for hydrolase acting on sulfur-nitrogen bonds.

    Examples include:
    - Sulfonamidases
    - N-sulfoglucosamine sulfohydrolase
    """

    EC_NUMBER_PREFIX = "3.10.-.-"  # EC prefix for S-N bond hydrolases

    @abstractmethod
    def _implementation_only(self) -> None:
        """Mark this legacy implementation class as abstract."""
    
    # Known S-N bond containing compounds
    SN_BOND_COMPOUNDS = {
        "CHEBI:57868",  # N-sulfo-D-glucosamine
        "CHEBI:41522",  # cyclohexylsulfamate (alternative ID)
        "CHEBI:57592",  # cyclohexylsulfamate (RHEA uses this ID)
        # Add more known S-N compounds as discovered
    }
    
    # Compounds that are sulfate esters (O-S bonds, not S-N)
    SULFATE_ESTER_COMPOUNDS = {
        "CHEBI:16189",  # sulfate (product of both S-N and O-S hydrolysis)
        # Note: sulfate alone doesn't indicate S-N vs O-S
    }
    
    # Name patterns that indicate S-N bonds
    SN_NAME_PATTERNS = [
        "n-sulfo",
        "n-sulfon",
        "sulfamate",
        "sulfonamide",
        "sulfamide",
        "sulfon-amid",
        "amino.*sulf",
        "sulf.*amino",
    ]
    
    # Name patterns that indicate O-S bonds (not S-N)
    OS_NAME_PATTERNS = [
        "sulfate ester",
        "o-sulfate",
        "o-sulfo",
        "-yl sulfate",
        "sulfate-4",
        "sulfate-6",
        "4-sulfate",
        "6-sulfate",
    ]

    # Define patterns using DSL
    # S-N bond hydrolysis patterns - more specific
    PATTERNS: list[Reaction] = [
        # N-sulfo compound + H2O → amine + sulfate
        # e.g., N-sulfo-D-glucosamine + H2O = D-glucosamine + sulfate
        var("n_sulfo_compound") + p(water)
        >> var("amine") + var("sulfate"),
        
        # Sulfamate + H2O → amine + sulfate  
        # e.g., cyclohexylsulfamate + H2O = cyclohexylamine + sulfate
        var("sulfamate") + p(water)
        >> var("amine") + var("sulfate"),
        
        # With H+ production: substrate + H2O → products + H+
        var("substrate") + p(water)
        >> var("product1") + var("product2") + h_plus,
        
        # General S-N hydrolysis with optional H+
        var("substrate_sn") + p(water)
        >> var("product_s") + var("product_n") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if a reaction is hydrolysis of S-N bonds using pattern matching.

        Strategy:
        1. Check for known S-N compounds by ChEBI ID
        2. Check compound names for S-N patterns  
        3. Exclude O-S bond reactions (sulfate esters)
        4. Basic hydrolase requirements (water consumption)
        5. Verify fragmentation pattern
        """
        # First check if it's even a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase reaction: {parent_result.explanation}",
            )
        
        # Check for known S-N bond compounds by ChEBI ID FIRST
        has_sn_compound = False
        sn_compound_name = None

        for participant in reaction.left_participants:
            if participant.chebi_id in self.SN_BOND_COMPOUNDS:
                has_sn_compound = True
                sn_compound_name = participant.name or participant.chebi_id
                break

        # Check by SMARTS-based sulfamate detection (structural)
        # But only if S-N bond is actually broken (not preserved in products)
        if not has_sn_compound:
            for participant in reaction.left_participants:
                if participant.is_sulfamate():
                    # Check if S-N bond is preserved in any product
                    # If a product also has sulfamate, it's O-S hydrolysis not S-N
                    sn_preserved_in_products = any(
                        prod.is_sulfamate() for prod in reaction.right_participants
                    )
                    if not sn_preserved_in_products:
                        has_sn_compound = True
                        sn_compound_name = participant.name or participant.chebi_id or "sulfamate (SMARTS)"
                        break

        # Check by name patterns as fallback (for molecules without SMILES)
        if not has_sn_compound:
            import re
            for participant in reaction.left_participants:
                if participant.name:
                    name_lower = participant.name.lower()
                    # Check for S-N bond patterns
                    for pattern in self.SN_NAME_PATTERNS:
                        if re.search(pattern, name_lower):
                            has_sn_compound = True
                            sn_compound_name = participant.name
                            break

                    if has_sn_compound:
                        break
        
        # Initialize variable that will be used later
        has_sulfate_product = False
        
        # If we have clear S-N compound, that's strong evidence
        if has_sn_compound:
            # Check for products that suggest S-N hydrolysis (sulfate)
            has_sulfate_product = any(
                participant.chebi_id == "CHEBI:16189" or 
                (participant.name and "sulfate" in participant.name.lower())
                for participant in reaction.right_participants
            )

        # Check for O-S patterns (sulfate ester) - these are NOT S-N bonds
        # Must check BEFORE fallback inference to avoid false positives
        import re
        has_os_pattern = False
        for participant in reaction.left_participants:
            if participant.name:
                name_lower = participant.name.lower()
                # Check if it matches O-S patterns (sulfate ester)
                for pattern in self.OS_NAME_PATTERNS:
                    if re.search(pattern, name_lower):
                        has_os_pattern = True
                        break
                if has_os_pattern:
                    break

        # If we have O-S pattern and no explicit S-N compound, it's sulfate ester hydrolysis
        if has_os_pattern and not has_sn_compound:
            return ClassificationResult(
                is_member=False,
                explanation="O-S bond hydrolysis (sulfate ester), not S-N bond"
            )

        # If no S-N compound found, check for S and N as fallback
        if not has_sn_compound:
            diff = ReactionDiff(reaction)
            has_sulfur = "S" in diff.reactant_elements
            has_nitrogen = "N" in diff.reactant_elements
            
            # If we have both S and N, this could be S-N hydrolysis
            if has_sulfur and has_nitrogen:
                # Check for sulfate product (common in S-N hydrolysis)
                has_sulfate_product = any(
                    participant.chebi_id == "CHEBI:16189" or 
                    (participant.name and "sulfate" in participant.name.lower())
                    for participant in reaction.right_participants
                )
                
                # Check for amine/ammonia product (from N part)
                has_nitrogen_product = any(
                    participant.name and ("amine" in participant.name.lower() or "amino" in participant.name.lower() or 
                               "ammonium" in participant.name.lower() or "ammonia" in participant.name.lower())
                    for participant in reaction.right_participants
                )
                
                # Require BOTH sulfate product AND nitrogen-containing substrate
                # to infer S-N bond hydrolysis (was too permissive with OR)
                if has_sulfate_product and has_nitrogen_product:
                    # Likely S-N hydrolysis even without explicit compound name
                    has_sn_compound = True
                    sn_compound_name = "S-N bond compound (inferred)"

            if not (has_sulfur and has_nitrogen):
                return ClassificationResult(
                    is_member=False, 
                    explanation="No S-N bond compound detected (missing S or N)"
                )

        # At this point we have S-N compound - check products
        if has_sn_compound:
            # Check for amine product (expected from S-N cleavage)
            has_amine_product = any(
                participant.name and ("amine" in participant.name.lower() or "amino" in participant.name.lower())
                for participant in reaction.right_participants
            )

            # Try pattern matching
            match = match_patterns(reaction, self.PATTERNS, strict=False)

            # Build explanation
            explanation = f"S-N bond hydrolysis: {sn_compound_name}"
            
            if has_sulfate_product:
                explanation += " → sulfate"
            if has_amine_product:
                explanation += " + amine"
            
            # Check for acid conditions
            if self._check_acid_conditions(reaction):
                explanation += ", acid conditions"
            
            if match and match.matched:
                explanation += " [pattern matched]"
            
            return ClassificationResult(is_member=True, explanation=explanation)

        # If we still don't have evidence, be conservative
        return ClassificationResult(
            is_member=False,
            explanation="Insufficient evidence for S-N bond hydrolysis",
        )

    def get_ec_number(self) -> str:
        """Return EC number prefix for S-N bond hydrolases."""
        return self.EC_NUMBER_PREFIX

    @staticmethod
    def _check_acid_conditions(reaction: Reaction) -> bool:
        """Check if reaction occurs under acid conditions.

        This is a simplified check - in reality would look for:
        - Protonated species
        - pH indicators
        - Acid catalysts

        For now, assume neutral/acid conditions if not explicitly basic
        (since the GO term specifies "acid" S-N bonds)
        """
        # SMILES-based acid detection would need RDKit - for now assume True
        return True
