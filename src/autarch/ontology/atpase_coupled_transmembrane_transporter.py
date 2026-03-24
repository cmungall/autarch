"""ATPase-coupled translocase reaction classification.

A more specific translocase that specifically requires ATP hydrolysis
as the energy source for active transport.

GO Hierarchy:
- GO:0015399 primary active transmembrane transporter activity (parent)
  - GO:0042626 ATPase-coupled transmembrane transporter activity (this class)
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_ADP,
)


class ATPaseCoupledTransmembraneTransporter(PrimaryActiveTransmembraneTransporter):
    """ATPase-coupled transmembrane transporter activity.

    A more specific translocase that specifically requires ATP hydrolysis
    as the energy source for active transport.

    This is a subclass of PrimaryActiveTransmembraneTransporter (GO:0015399) for reactions
    specifically coupled to ATP hydrolysis.

    Examples:
    - Ca2+-ATPase: Ca2+(in) + ATP + H2O → Ca2+(out) + ADP + Pi + H+
    - Na+/K+-ATPase: 3Na+(in) + 2K+(out) + ATP → 3Na+(out) + 2K+(in) + ADP + Pi
    - ABC transporters: substrate(out) + ATP → substrate(in) + ADP + Pi
    """

    GO_ID: ClassVar[Optional[str]] = "GO:0042626"  # ATPase-coupled transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = None  # No specific EC prefix, use parent's

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is specifically an ATPase-coupled translocase.

        Strategy:
        1. Must pass parent PrimaryActiveTransmembraneTransporter checks
        2. Must specifically use ATP (not other energy sources)
        """
        # First check parent class
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        # Must specifically use ATP
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )

        has_adp = any(
            p.chebi_id == CHEBI_ADP
            for p in reaction.right_participants
        )

        if not (has_atp and has_adp):
            return ClassificationResult(
                is_member=False,
                explanation="Not ATP-coupled - uses other energy source"
            )

        # Get transported molecules for explanation
        transported = reaction.get_transported_molecules()
        transport_info = []
        for mol_id, from_loc, to_loc in transported:
            transport_info.append(f"{mol_id}: {from_loc}→{to_loc}")

        return ClassificationResult(
            is_member=True,
            explanation=f"ATPase-coupled translocase: {', '.join(transport_info)}"
        )
