"""Peptidyl-prolyl Isomerase reaction classification using pattern DSL.

Peptidyl-prolyl isomerases (PPIases) are ultra-specific isomerases that catalyze
the cis-trans isomerization of proline residues in peptide chains. This is a
unique protein folding assistance mechanism.

Key biochemical signature: peptidyl-L-proline (cis) ⇌ peptidyl-L-proline (trans)

This tests our ability to classify highly specific protein folding catalysts
with no cofactor requirements - pure conformational catalysis.

History
-------

## 2025-12-22

Fixed critical bug: classifier was not calling parent Isomerase.check_membership_impl(),
causing it to match any 1-2 reactant/product reaction without cofactors. Now properly
inherits from Isomerase and adds proline-specific detection.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_NADH,
    CHEBI_NADPH,
    water,
)


class PeptidylProlylCisTransIsomerase(Isomerase):
    """peptidyl-prolyl cis-trans isomerase

    Examples:
    - Cyclophilin: immunosuppressant-binding PPIase
    - FKBP (FK506-binding protein): tacrolimus-sensitive PPIase
    - Pin1: phospho-specific PPIase for Ser/Thr-Pro motifs

    This is ULTRA-SPECIFIC:
    - Only ~10-15 reactions in databases
    - Requires proline residues in peptides
    - No cofactors needed (pure conformational catalysis)
    - Essential for protein folding
    """

    GO_ID = "GO:0003755"  # peptidyl-prolyl cis-trans isomerase activity
    EC_NUMBER_PREFIX = "5.2.1.8"  # Peptidylprolyl isomerase

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple isomerization: peptidyl-proline (cis) ⇌ peptidyl-proline (trans)
        var("cis_peptidyl_proline") >> var("trans_peptidyl_proline"),
        # With water involvement (rare)
        var("cis_peptidyl_proline") + optional(water) >> var("trans_peptidyl_proline") + optional(water),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is peptidyl-prolyl isomerase using PROLINE-SPECIFIC criteria.

        PEPTIDYL-PROLYL ISOMERASE SIGNATURE:
        1. Must be an isomerase first (parent class)
        2. Should involve peptidyl-proline substrates
        3. Should be reversible isomerization (1→1 typically)
        4. Should NOT require cofactors (pure conformational)

        This tests recognition of conformational catalysis!
        """
        # First check if it's an isomerase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an isomerase: {parent_result.explanation}",
            )

        # Should NOT use major cofactors (ATP, NAD, etc.)
        has_major_cofactors = any(
            p.chebi_id in [CHEBI_ATP, CHEBI_NADPH, CHEBI_NADH]
            for p in (reaction.left_participants + reaction.right_participants)
        )

        if has_major_cofactors:
            return ClassificationResult(
                is_member=False,
                explanation="Uses major cofactors - PPIases are cofactor-independent"
            )

        # Known peptidyl-prolyl isomerase ChEBI IDs
        ppiase_chebis = {
            "CHEBI:83155",  # peptidyl-L-proline (omega = 0)
            "CHEBI:83154",  # peptidyl-L-proline (omega = 180)
        }

        # Check for proline-related substrates
        has_proline_substrate = any(
            p.chebi_id in ppiase_chebis
            for p in reaction.left_participants
        )

        has_proline_product = any(
            p.chebi_id in ppiase_chebis
            for p in reaction.right_participants
        )

        label_lower = reaction.label.lower() if reaction.label else ""
        has_proline_label = (
            "peptidylproline" in label_lower or
            "peptidyl-proline" in label_lower
        )

        if not (has_proline_substrate or has_proline_product or has_proline_label):
            return ClassificationResult(
                is_member=False,
                explanation="No peptidyl-proline substrates detected"
            )

        return ClassificationResult(
            is_member=True,
            explanation="Peptidyl-prolyl cis-trans isomerase: proline cis-trans isomerization"
        )
