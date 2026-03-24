"""Translocase linked to hydrolysis of a nucleoside triphosphate.

EC 7.2 translocases are the most common type, using ATP or GTP hydrolysis
to power active transport across membranes (ABC transporters, P-type ATPases).
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_ADP,
    CHEBI_GTP,
    CHEBI_GDP,
)


class TranslocaseLinkedToHydrolysis(PrimaryActiveTransmembraneTransporter):
    """ATPase-coupled transmembrane transporter activity

    Enables the transfer of a solute across a membrane, powered by
    hydrolysis of a nucleoside triphosphate (ATP or GTP).

    Examples:
    - ABC transporters (ATP-binding cassette)
    - P-type ATPases (Na+/K+-ATPase, Ca2+-ATPase)
    - V-type ATPases (vacuolar proton pumps)
    - F-type ATPases (in reverse, ATP-driven proton transport)
    """

    GO_ID: ClassVar[Optional[str]] = "GO:0042626"  # ATPase-coupled transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.2.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an NTP hydrolysis-coupled translocase.

        Strategy:
        1. Must be a transport reaction (from parent)
        2. Must have NTP (ATP/GTP) as reactant AND NDP (ADP/GDP) as product
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

        # Check for NTP hydrolysis (ATP -> ADP or GTP -> GDP)
        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        has_atp_hydrolysis = CHEBI_ATP in left_ids and CHEBI_ADP in right_ids
        has_gtp_hydrolysis = CHEBI_GTP in left_ids and CHEBI_GDP in right_ids

        if not (has_atp_hydrolysis or has_gtp_hydrolysis):
            return ClassificationResult(
                is_member=False,
                explanation="No NTP hydrolysis (ATP->ADP or GTP->GDP) - not hydrolysis-coupled translocase",
            )

        ntp_type = "ATP" if has_atp_hydrolysis else "GTP"
        transport_info = [f"{mol_id}: {f}→{t}" for mol_id, f, t in transported]
        return ClassificationResult(
            is_member=True,
            explanation=f"Translocase ({ntp_type} hydrolysis-coupled): {', '.join(transport_info)}",
        )
