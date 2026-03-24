"""Oxidoreductase forming X-Y bond reaction classification.

EC 1.21: Oxidoreductases that catalyze X-H + Y-H -> X-Y reactions.
These enzymes oxidize two substrates by forming a new bond between them,
with the removal of hydrogen. Examples include iodide peroxidase
(thyroid peroxidase), sulfhydryl oxidase (disulfide bond formation),
and thiol oxidase.

Key features:
- Bond formation between two donor substrates
- Iodide peroxidase (iodination of thyroglobulin)
- Sulfhydryl oxidase (thiol -> disulfide)
- Often uses O2 or H2O2 as acceptor
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import CHEBI_O2, CHEBI_H2O2


class OxidoreductaseFormingXYBond(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which hydrogen or electrons are transferred from each of two donors, and the two donors become joined."""

    GO_ID = "GO:0046992"  # oxidoreductase activity, forming X-Y bonds
    EC_NUMBER_PREFIX = "1.21.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase that forms X-Y bonds.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Label-based detection for known EC 1.21 enzymes
        3. Look for iodide peroxidase, sulfhydryl oxidase patterns
        """
        # First check parent
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        # Label-based detection for known EC 1.21 reactions
        xy_bond_labels = [
            "iodide peroxidase",
            "sulfhydryl oxidase",
            "thiol oxidase",
            "iodotyrosine",
            "diiodotyrosine",
        ]
        has_label_match = any(pat in label_lower for pat in xy_bond_labels)

        if has_label_match:
            return ClassificationResult(
                is_member=True,
                explanation=f"X-Y bond-forming oxidoreductase: label match in '{label_lower}'",
            )

        # Check for iodide involvement (iodide peroxidase)
        has_iodide = any(
            p.name and "iodide" in p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
        ) or "iodide" in label_lower

        if has_iodide:
            # Must also have H2O2 or O2 (peroxidase mechanism)
            has_oxidant = any(
                p.chebi_id in {CHEBI_O2, CHEBI_H2O2}
                for p in reaction.left_participants + reaction.right_participants
            )
            if has_oxidant:
                return ClassificationResult(
                    is_member=True,
                    explanation="X-Y bond-forming oxidoreductase: iodide with oxidant",
                )

        return ClassificationResult(
            is_member=False,
            explanation="No X-Y bond-forming oxidoreductase pattern detected",
        )
