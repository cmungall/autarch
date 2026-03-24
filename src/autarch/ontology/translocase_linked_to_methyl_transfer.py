"""Translocase linked to methyl transfer.

EC 7.4 translocases use methyl group transfer energy to power active transport.
These are extremely rare enzymes.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.molecules import CHEBI_SAM, CHEBI_SAH


class TranslocaseLinkedToMethylTransfer(PrimaryActiveTransmembraneTransporter):
    """primary active transmembrane transporter activity

    Translocase powered by methyl group transfer.
    The energy released by transferring a methyl group
    (e.g., from S-adenosylmethionine) drives transport across the membrane.

    Examples:
    - Na+-transporting methyltransferase (methanogenesis pathway)
    """

    GO_ID = "GO:0015399"  # primary active transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.4.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a methyl transfer-coupled translocase.

        Strategy:
        1. Must be a transport reaction
        2. Must show methyl transfer (SAM -> SAH or similar)
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

        # Check for methyl transfer (SAM consumed, SAH produced)
        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        has_methyl_transfer = CHEBI_SAM in left_ids and CHEBI_SAH in right_ids

        # Also check labels for methyl transfer indicators
        all_names = [
            p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
            if p.name
        ]
        has_methyl_label = any("methyl" in name for name in all_names)

        if not (has_methyl_transfer or has_methyl_label):
            return ClassificationResult(
                is_member=False,
                explanation="No methyl transfer detected - not methyl transfer-coupled translocase",
            )

        transport_info = [f"{mol_id}: {f}→{t}" for mol_id, f, t in transported]
        return ClassificationResult(
            is_member=True,
            explanation=f"Translocase (methyl transfer-coupled): {', '.join(transport_info)}",
        )
