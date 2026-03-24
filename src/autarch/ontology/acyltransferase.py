"""Acyltransferase reaction classification using pattern DSL.

Acyltransferases transfer acyl groups from donors to acceptors.
EC 2.3.x.x classification.

History
-------

## 2025-12-21

Refactored to use SMARTS-based thioester detection (Moiety.THIOESTER) for CoA
derivative recognition instead of name-based matching. Added Participant.is_thioester()
check alongside ChEBI ID matching for acyl donor detection.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase
from autarch.pattern_dsl import var, match_patterns
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
)

# Electron transfer flavoprotein (ETF) for acyl-CoA dehydrogenase exclusion
CHEBI_ETF_OX = "CHEBI:57692"   # oxidized [electron-transfer flavoprotein]
CHEBI_ETF_RED = "CHEBI:57716"  # reduced [electron-transfer flavoprotein]


class Acyltransferase(Transferase):
    """acyltransferase

    Examples:
    - Acetyl-CoA acetyltransferase: 2 acetyl-CoA → acetoacetyl-CoA + CoA
    - Choline acetyltransferase: acetyl-CoA + choline → acetylcholine + CoA
    - Chloramphenicol acetyltransferase: acetyl-CoA + chloramphenicol → N-acetylchloramphenicol + CoA
    - Fatty acid synthase: acetyl-CoA + malonyl-CoA → fatty acid intermediates
    """

    GO_ID = "GO:0016746"  # acyltransferase activity
    EC_NUMBER_PREFIX = "2.3.-.-"  # Acyltransferases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple acyl transfer: acyl-CoA + acceptor → product + CoA
        var("acyl_coa") + var("acceptor") >> var("product") + var("coa"),
        # Acetyl transfer: acetyl-CoA + acceptor → acetyl-acceptor + CoA
        var("acetyl_coa") + var("acceptor") >> var("acetyl_product") + var("coa"),
        # Condensation: acyl-CoA + acyl-CoA → longer-acyl-CoA + CoA
        var("acyl_coa1") + var("acyl_coa2") >> var("condensed_product") + var("coa"),
        # With cofactor: acyl-donor + acceptor + cofactor → product + modified-donor + modified-cofactor
        var("donor") + var("acceptor") + var("cofactor") >> var("product") + var("modified_donor") + var("modified_cofactor"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an acyltransferase using pattern matching.

        Strategy:
        1. Must be a transferase (parent class)
        2. Must involve acyl group transfer (usually CoA derivatives)
        3. Look for CoA-dependent acyl transfer patterns
        4. Exclude non-acyl transfers
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's a transferase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a transferase: {parent_result.explanation}",
            )

        # Look for acyl group donors (typically CoA derivatives)
        # Use both ChEBI ID matching and SMARTS-based thioester detection
        acyl_donors = {
            CHEBI_COA,  # acetyl-CoA
            "CHEBI:15351",  # malonyl-CoA
            "CHEBI:57804",  # palmitoyl-CoA
            "CHEBI:15525",  # propionyl-CoA
            "CHEBI:57262",  # butyryl-CoA
        }

        # Check for acyl donors via ChEBI ID or SMARTS-based thioester detection
        has_acyl_donor = any(
            p.chebi_id in acyl_donors or p.is_thioester()
            for p in reaction.left_participants
        )

        # Check for CoA as product (released when acyl group is transferred)
        has_coa_product = any(
            p.chebi_id == CHEBI_COA
            for p in reaction.right_participants
        )

        # Must have either acyl donor or CoA release pattern
        if not (has_acyl_donor or has_coa_product):
            return ClassificationResult(
                is_member=False,
                explanation="No acyl group transfer pattern detected"
            )

        # Apply exclusions before classifying as acyltransferase

        # Exclude oxidative reactions (NAD+/NADP+ involvement)
        has_nad_system = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}for p in reaction.left_participants + reaction.right_participants
        )

        if has_nad_system:
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P) system present - likely oxidoreductase not acyltransferase"
            )

        # Exclude acyl-CoA dehydrogenases (use ETF as electron acceptor)
        # Pattern: acyl-CoA + oxidized ETF → enoyl-CoA + reduced ETF
        etf_patterns = ["electron-transfer flavoprotein", "etf", "[etf]"]
        has_etf = any(
            p.chebi_id in {CHEBI_ETF_OX, CHEBI_ETF_RED}
            for p in reaction.left_participants + reaction.right_participants
        ) or any(pat in label_lower for pat in etf_patterns)
        if has_etf:
            return ClassificationResult(
                is_member=False,
                explanation="ETF cofactor present - acyl-CoA dehydrogenase not acyltransferase"
            )

        # Exclude CoA-transferases (EC 2.8.3.x) - transfer CoA group, not acyl group
        # Pattern: acyl-CoA + acid → CoA-acid + carboxylate (same CoA, different acyl)
        # Key indicator: succinyl-CoA + substrate → product-CoA + succinate
        SUCCINYL_COA_CHEBIS = {"CHEBI:15380", "CHEBI:57292"}  # succinyl-CoA
        SUCCINATE_CHEBIS = {"CHEBI:15741", "CHEBI:30031"}  # succinate
        has_succinyl_coa = any(
            p.chebi_id in SUCCINYL_COA_CHEBIS
            for p in reaction.left_participants
        ) or "succinyl-coa" in substrate_str
        has_succinate = any(
            p.chebi_id in SUCCINATE_CHEBIS
            for p in reaction.right_participants
        ) or "succinate" in product_str
        # Also check for acetyl-CoA + acid → acyl-CoA + acetate pattern
        ACETATE_CHEBIS = {"CHEBI:30089", "CHEBI:40480"}  # acetate
        ACETYL_COA_CHEBIS = {"CHEBI:57288", "CHEBI:15351"}  # acetyl-CoA
        has_acetyl_coa_donor = any(
            p.chebi_id in ACETYL_COA_CHEBIS
            for p in reaction.left_participants
        ) or "acetyl-coa" in substrate_str
        has_acetate_product = any(
            p.chebi_id in ACETATE_CHEBIS
            for p in reaction.right_participants
        ) or "acetate" in product_str
        # Check for non-CoA acyl-CoA product (indicates CoA transfer)
        # Use label pattern: "-coa" in product but not free "coa"
        has_other_acyl_coa_product = "-coa" in product_str and "coa +" not in product_str
        # CoA-transferase: CoA is transferred from one acyl to another
        if (has_succinyl_coa and has_succinate and has_other_acyl_coa_product):
            return ClassificationResult(
                is_member=False,
                explanation="CoA-transferase - transfers CoA group not acyl group"
            )
        if (has_acetyl_coa_donor and has_acetate_product and has_other_acyl_coa_product):
            return ClassificationResult(
                is_member=False,
                explanation="CoA-transferase - transfers CoA group not acyl group"
            )

        # Exclude simple hydrolysis (water + substrate → products + CoA)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )
        
        # Count non-water, non-CoA substrates
        non_cofactor_substrates = [
            p for p in reaction.left_participants
            if p.chebi_id not in {CHEBI_H2O, CHEBI_COA, CHEBI_H_PLUS} and  # H2O, CoA, H+
            not False ]
        
        if has_water and len(non_cofactor_substrates) == 1 and has_coa_product:
            return ClassificationResult(
                is_member=False,
                explanation="Simple hydrolysis - not acyl transfer"
            )

        # Exclude ligase reactions (ATP consumption with bond formation)
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        has_amp_or_adp = any(
            p.chebi_id in {CHEBI_AMP, CHEBI_ADP}
            for p in reaction.right_participants
        )
        
        if has_atp and has_amp_or_adp and len(non_cofactor_substrates) >= 2:
            return ClassificationResult(
                is_member=False,
                explanation="ATP-dependent ligase - not simple acyltransferase"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Determine specific type
        explanation = "Acyltransferase: acyl group transfer"
        
        # Check for specific acyl groups
        # Specific acyltransferase type detection would need ChEBI IDs

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
