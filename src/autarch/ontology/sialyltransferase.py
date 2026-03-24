"""sialyltransferase.

Catalysis of the transfer of a sialyl group to an acceptor.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_CMP, CHEBI_H2O, CHEBI_H_PLUS
from autarch.ontology.glycosyltransferase import Glycosyltransferase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE

CHEBI_CMP_NEUAC = "CHEBI:57812"


class Sialyltransferase(Glycosyltransferase):
    """sialyltransferase.

    Catalysis of the transfer of a sialyl group to an acceptor.
    """

    GO_ID = "GO:0008373"
    EC_NUMBER_PREFIX = "2.4.3.-"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    LEFT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = {CHEBI_H2O, CHEBI_H_PLUS}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for CMP-sialic-acid-dependent glycosyl transfer."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not sialyl transfer",
            )

        forward = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse.explanation} (reverse reaction orientation)",
            )

        return forward

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        has_sialyl_donor = any(
            participant.chebi_id == CHEBI_CMP_NEUAC for participant in left_participants
        )
        if not has_sialyl_donor:
            return ClassificationResult(
                is_member=False,
                explanation="Missing CMP-sialic-acid donor on the substrate side",
            )

        cmp_products = sum(
            max(participant.count, 1)
            for participant in right_participants
            if participant.chebi_id == CHEBI_CMP
        )
        if cmp_products == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Missing CMP leaving-group product expected for sialyl transfer",
            )

        acceptors = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
            and participant.chebi_id != CHEBI_CMP_NEUAC
        ]
        if not acceptors:
            return ClassificationResult(
                is_member=False,
                explanation="No acceptor substrate remained after removing the CMP-sialic-acid donor and solvent spectators",
            )

        transferred_products = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
            and participant.chebi_id != CHEBI_CMP
        ]
        if not transferred_products:
            return ClassificationResult(
                is_member=False,
                explanation="No glycosylated product remained after removing CMP and solvent spectators",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Sialyltransferase: a CMP-sialic-acid donor transfers a sialyl group to an acceptor with CMP release",
        )
