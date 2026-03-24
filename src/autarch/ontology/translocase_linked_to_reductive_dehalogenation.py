"""Translocase linked to reductive dehalogenation.

EC 7.6 translocases use reductive dehalogenation energy to power active transport.
These are extremely rare enzymes found in organohalide-respiring bacteria.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.primary_active_transmembrane_transporter import (
    PrimaryActiveTransmembraneTransporter,
)


# Common halide product IDs
CHEBI_CHLORIDE = "CHEBI:17996"
CHEBI_BROMIDE = "CHEBI:15858"
CHEBI_FLUORIDE = "CHEBI:17051"
CHEBI_IODIDE = "CHEBI:16382"

HALIDE_IDS = {CHEBI_CHLORIDE, CHEBI_BROMIDE, CHEBI_FLUORIDE, CHEBI_IODIDE}


class TranslocaseLinkedToReductiveDehalogenation(PrimaryActiveTransmembraneTransporter):
    """primary active transmembrane transporter activity

    Translocase powered by reductive dehalogenation.
    The energy released by removing a halogen substituent (with reduction)
    drives transport across the membrane.

    Examples:
    - Organohalide-respiring bacteria use reductive dehalogenation
      coupled to proton translocation across membranes
    """

    GO_ID = "GO:0015399"  # primary active transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.6.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a reductive dehalogenation-coupled translocase.

        Strategy:
        1. Must be a transport reaction
        2. Must produce halide ion (dehalogenation)
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

        # Check for halide as product (dehalogenation)
        has_halide_product = any(
            p.chebi_id in HALIDE_IDS
            for p in reaction.right_participants
        )

        # Also check labels for halide/dehalogenation indicators
        all_names = [
            p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
            if p.name
        ]
        has_halide_label = any(
            term in name
            for name in all_names
            for term in ("chloride", "bromide", "fluoride", "iodide", "dehalogen")
        )

        if not (has_halide_product or has_halide_label):
            return ClassificationResult(
                is_member=False,
                explanation="No halide produced - not reductive dehalogenation-coupled translocase",
            )

        transport_info = [f"{mol_id}: {f}→{t}" for mol_id, f, t in transported]
        return ClassificationResult(
            is_member=True,
            explanation=f"Translocase (reductive dehalogenation-coupled): {', '.join(transport_info)}",
        )
