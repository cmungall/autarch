"""Translocase linked to deamination.

EC 7.5 translocases use deamination energy to power active transport.
These are extremely rare enzymes.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.molecules import CHEBI_NH3, CHEBI_NH4


class TranslocaseLinkedToDeamination(PrimaryActiveTransmembraneTransporter):
    """primary active transmembrane transporter activity

    Translocase powered by deamination of a substrate.
    The energy released by removing an amino group as ammonia
    drives transport across the membrane.

    Examples:
    - Hypothetical deamination-coupled Na+ pump
    """

    GO_ID = "GO:0015399"  # primary active transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.5.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a deamination-coupled translocase.

        Strategy:
        1. Must be a transport reaction
        2. Must produce ammonia/ammonium (deamination)
        3. Does NOT require ATP
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

        # Check for ammonia/ammonium as product (deamination)
        has_ammonia_product = any(
            p.chebi_id in {CHEBI_NH3, CHEBI_NH4}
            for p in reaction.right_participants
        )

        if not has_ammonia_product:
            return ClassificationResult(
                is_member=False,
                explanation="No ammonia produced - not deamination-coupled translocase",
            )

        transport_info = [f"{mol_id}: {f}→{t}" for mol_id, f, t in transported]
        return ClassificationResult(
            is_member=True,
            explanation=f"Translocase (deamination-coupled): {', '.join(transport_info)}",
        )
