"""Translocase linked to decarboxylation.

EC 7.3 translocases use decarboxylation energy to power active transport.
These are rare enzymes, primarily sodium-pumping decarboxylases.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.molecules import CHEBI_CO2


class TranslocaseLinkedToDecarboxylation(PrimaryActiveTransmembraneTransporter):
    """primary active transmembrane transporter activity

    Translocase powered by decarboxylation of a substrate.
    The energy released by removing a carboxyl group as CO2
    drives transport across the membrane.

    Examples:
    - Oxaloacetate decarboxylase Na+ pump
    - Methylmalonyl-CoA decarboxylase Na+ pump
    """

    GO_ID = "GO:0015451"  # decarboxylation-driven active transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.3.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a decarboxylation-coupled translocase.

        Strategy:
        1. Must be a transport reaction
        2. Must produce CO2 (decarboxylation)
        3. Does NOT require ATP (unlike parent class)
        """
        # Must be a transport reaction
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="No molecule transported between locations",
            )

        transported = reaction.get_transported_molecules()
        if not transported:
            return ClassificationResult(
                is_member=False,
                explanation="No transported molecules detected",
            )

        # Check for CO2 as product (decarboxylation)
        has_co2_product = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.right_participants
        )

        if not has_co2_product:
            return ClassificationResult(
                is_member=False,
                explanation="No CO2 produced - not decarboxylation-coupled translocase",
            )

        transport_info = [f"{mol_id}: {f}→{t}" for mol_id, f, t in transported]
        return ClassificationResult(
            is_member=True,
            explanation=f"Translocase (decarboxylation-coupled): {', '.join(transport_info)}",
        )
