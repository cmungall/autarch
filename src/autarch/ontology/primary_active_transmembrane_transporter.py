"""Translocase reaction classification.

Translocases catalyze active transport of ions or molecules across membranes.
EC 7.x.x.x in the enzyme classification system (added in 2018).

GO Hierarchy:
- GO:0015399 primary active transmembrane transporter activity (= EC:7.-.-.-)
  - GO:0042626 ATPase-coupled transmembrane transporter activity (more specific)

Detection is based on:
- Location annotations showing same molecule on both sides (in/out)
- Primary energy source coupling (ATP, redox, or photon energy)
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_ADP,
    CHEBI_NADH,
    CHEBI_NAD_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
)


class PrimaryActiveTransmembraneTransporter(ReactionClass):
    """primary active transmembrane transporter

    Enables the transfer of a solute from one side of a membrane to the other,
    up the solute's concentration gradient. Transport is powered by a primary
    energy source (chemical such as ATP hydrolysis, redox energy, or photon energy).

    This is the generic EC 7.x.x.x translocase class.

    Common examples:
    - Ion pumps driven by ATP, redox, or light
    - ABC transporters
    - P-type ATPases

    Note: EC 7.x.x.x is a newer classification (2018) for translocases.
    Many are still classified under EC 3.6.3.x (ATP-driven transport).
    """

    GO_ID: ClassVar[Optional[str]] = "GO:0015399"  # primary active transmembrane transporter activity
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "7.-.-.-"  # Translocases
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a primary active translocase.

        Strategy:
        1. Must be a transport reaction (same molecule, different locations)
        2. Must have primary energy coupling (ATP, redox, or similar)
        """
        # Must be a transport reaction (same molecule, different locations)
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="No molecule transported between locations"
            )

        # Get transported molecules
        transported = reaction.get_transported_molecules()
        if not transported:
            return ClassificationResult(
                is_member=False,
                explanation="No transported molecules detected"
            )

        # Check for primary energy sources
        # ATP coupling
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )

        # Redox coupling (NAD(P)H consumption)
        has_redox = any(
            p.chebi_id in {CHEBI_NADH, CHEBI_NADPH}
            for p in reaction.left_participants
        )

        # Check for products indicating energy consumption
        has_adp = any(
            p.chebi_id == CHEBI_ADP
            for p in reaction.right_participants
        )

        has_oxidized_cofactor = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
            for p in reaction.right_participants
        )

        # Must have primary energy coupling
        has_energy_coupling = (has_atp and has_adp) or (has_redox and has_oxidized_cofactor) or has_atp

        if not has_energy_coupling:
            return ClassificationResult(
                is_member=False,
                explanation="No primary energy coupling - passive transport, not translocase"
            )

        # Build explanation with transport details
        transport_info = []
        for mol_id, from_loc, to_loc in transported:
            transport_info.append(f"{mol_id}: {from_loc}→{to_loc}")

        energy_type = "ATP" if has_atp else "redox" if has_redox else "primary energy"
        return ClassificationResult(
            is_member=True,
            explanation=f"Primary active translocase ({energy_type}-coupled): {', '.join(transport_info)}"
        )
