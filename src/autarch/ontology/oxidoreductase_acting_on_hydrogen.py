"""Oxidoreductase acting on hydrogen as donor.

EC 1.12: Hydrogenases that use molecular hydrogen (H2) as electron donor.

These enzymes catalyze the oxidation of H2:
- [NiFe]-hydrogenase: H2 + acceptor → 2H+ + reduced acceptor
- [FeFe]-hydrogenase: H2 + ferredoxin(ox) → 2H+ + ferredoxin(red)
- H2:NAD+ oxidoreductase: H2 + NAD+ → H+ + NADH
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Molecular hydrogen
CHEBI_H2 = "CHEBI:18276"  # molecular hydrogen (dihydrogen)


class OxidoreductaseActingOnHydrogen(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which hydrogen acts as an electron donor.

    EC 1.12.x.x includes:
    - [NiFe]-hydrogenase: H2 + acceptor = 2 H+ + reduced acceptor
    - [FeFe]-hydrogenase: H2 + ferredoxin(ox) = 2 H+ + ferredoxin(red)
    - H2:NAD+ oxidoreductase: H2 + NAD+ = H+ + NADH

    Examples:
    - H2 + NAD+ = H+ + NADH
    - H2 + 2 ferricytochrome c3 = 2 H+ + 2 ferrocytochrome c3
    - H2 + acceptor = reduced acceptor
    """

    GO_ID = "GO:0016695"  # oxidoreductase activity, acting on hydrogen as donor
    EC_NUMBER_PREFIX = "1.12.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a hydrogenase (H2 as donor).

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Must have molecular hydrogen (H2) as a reactant
        This is a very specific class: H2 is the defining substrate.
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Check for molecular hydrogen as reactant (donor)
        has_h2_reactant = any(
            p.chebi_id == CHEBI_H2
            for p in reaction.left_participants
        )

        if has_h2_reactant:
            return ClassificationResult(
                is_member=True,
                explanation="Hydrogenase: molecular hydrogen (H2) as electron donor"
            )

        # Also check for H2 as product (reverse hydrogenase / H2-evolving)
        has_h2_product = any(
            p.chebi_id == CHEBI_H2
            for p in reaction.right_participants
        )

        if has_h2_product:
            return ClassificationResult(
                is_member=True,
                explanation="Hydrogenase: molecular hydrogen (H2) as product (reverse direction)"
            )

        # Check label for "dihydrogen" as a standalone word (the chemical name for H2)
        # Must avoid matching compound names like "dihydrogenistein"
        label_lower = reaction.label.lower() if reaction.label else ""
        if re.search(r'\bdihydrogen\b', label_lower):
            return ClassificationResult(
                is_member=True,
                explanation="Hydrogenase: label contains dihydrogen"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No molecular hydrogen (H2) detected as reactant or product"
        )
