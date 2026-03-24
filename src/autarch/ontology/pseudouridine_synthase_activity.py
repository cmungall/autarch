"""pseudouridine synthase activity.

Catalysis of the reaction: a uridine in RNA = a pseudouridine in RNA. Conversion of uridine in an RNA molecule to pseudouridine by rotation of the C1'-N-1 glycosidic bond of uridine in RNA to a C1'-C5.
"""

from __future__ import annotations

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_URIDINE_IN_RNA = "CHEBI:65315"
CHEBI_PSEUDOURIDINE_IN_RNA = "CHEBI:65314"


class PseudouridineSynthaseActivity(ReactionClass):
    """pseudouridine synthase activity.

    Catalysis of the reaction: a uridine in RNA = a pseudouridine in RNA. Conversion of uridine in an RNA molecule to pseudouridine by rotation of the C1'-N-1 glycosidic bond of uridine in RNA to a C1'-C5.
    """

    GO_ID = "GO:0009982"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = {
            participant.chebi_id
            for participant in reaction.left_participants
            if participant.chebi_id
        }
        right_ids = {
            participant.chebi_id
            for participant in reaction.right_participants
            if participant.chebi_id
        }

        if left_ids == {CHEBI_URIDINE_IN_RNA} and right_ids == {CHEBI_PSEUDOURIDINE_IN_RNA}:
            return ClassificationResult(
                is_member=True,
                explanation="Pseudouridine synthase activity: uridine in RNA isomerized to pseudouridine in RNA",
            )
        if left_ids == {CHEBI_PSEUDOURIDINE_IN_RNA} and right_ids == {CHEBI_URIDINE_IN_RNA}:
            return ClassificationResult(
                is_member=True,
                explanation="Pseudouridine synthase activity: uridine in RNA isomerized to pseudouridine in RNA (reverse reaction orientation)",
            )

        return ClassificationResult(
            is_member=False,
            explanation="Requires the supported uridine-in-RNA to pseudouridine-in-RNA isomerization pair",
        )
