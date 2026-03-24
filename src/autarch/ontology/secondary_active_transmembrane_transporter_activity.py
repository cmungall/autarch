"""secondary active transmembrane transporter activity."""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ADP, CHEBI_ATP, CHEBI_GDP, CHEBI_GTP, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.transport_utils import TRANSPORT_SPECTATOR_IDS, transported_pairs


class SecondaryActiveTransmembraneTransporterActivity(ReactionClass):
    """secondary active transmembrane transporter activity."""

    GO_ID = "GO:0015291"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    PRIMARY_ENERGY_IDS = {
        CHEBI_ATP,
        CHEBI_ADP,
        CHEBI_GTP,
        CHEBI_GDP,
        CHEBI_NADH,
        CHEBI_NADPH,
        CHEBI_NAD_PLUS,
        CHEBI_NADP_PLUS,
    }
    COUPLING_SPECTATOR_IDS = TRANSPORT_SPECTATOR_IDS - {CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate between locations",
            )

        transported = [
            (left_participant, right_participant)
            for left_participant, right_participant in transported_pairs(reaction)
            if left_participant.chebi_id not in self.COUPLING_SPECTATOR_IDS
        ]
        if len(transported) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Secondary active transport requires coupled substrate transport",
            )

        if any(
            participant.chebi_id in self.PRIMARY_ENERGY_IDS
            for participant in reaction.all_participants()
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Primary-energy coupling detected",
            )

        directions = {
            (left_participant.location, right_participant.location)
            for left_participant, right_participant in transported
        }
        mode = "symport" if len(directions) == 1 else "antiport"
        return ClassificationResult(
            is_member=True,
            explanation=f"Secondary active transmembrane transporter activity: coupled {mode}",
        )
