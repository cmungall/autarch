"""hydroxymethyl-, formyl- and related transferase activity.

Catalysis of the transfer of a hydroxymethyl- or formyl group from one
compound (donor) to another (acceptor).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.transferase import Transferase


class HydroxymethylFormylAndRelatedTransferase(Transferase):
    """hydroxymethyl-, formyl- and related transferase activity.

    Catalysis of the transfer of a hydroxymethyl- or formyl group from one
    compound (donor) to another (acceptor).
    """

    GO_ID = "GO:0016742"
    EC_NUMBER_PREFIX = "2.1.2.-"
    ONE_CARBON_FOLATE_DONORS = {
        "CHEBI:15636",
        "CHEBI:195366",
        "CHEBI:57456",
        "CHEBI:57457",
    }
    TETRAHYDROFOLATE_FORMS = {"CHEBI:57453"}

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for one-carbon transfer using tetrahydrofolate-derived carriers."""
        participant_chebis = {
            participant.chebi_id
            for participant in reaction.all_participants()
            if participant.chebi_id
        }
        if CHEBI_SAM in participant_chebis or CHEBI_SAH in participant_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="SAM/SAH methyl-transfer chemistry is not a folate one-carbon transferase",
            )
        if CHEBI_ATP in participant_chebis:
            return ClassificationResult(
                is_member=False,
                explanation="ATP-dependent chemistry is not hydroxymethyl/formyl transferase activity",
            )
        if participant_chebis & {CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH}:
            return ClassificationResult(
                is_member=False,
                explanation="NAD(P)-dependent redox chemistry is not hydroxymethyl/formyl transferase activity",
            )

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

        has_forward_folate_cycle = bool(left_ids & self.ONE_CARBON_FOLATE_DONORS) and bool(
            right_ids & self.TETRAHYDROFOLATE_FORMS
        )
        has_reverse_folate_cycle = bool(left_ids & self.TETRAHYDROFOLATE_FORMS) and bool(
            right_ids & self.ONE_CARBON_FOLATE_DONORS
        )
        if not (has_forward_folate_cycle or has_reverse_folate_cycle):
            return ClassificationResult(
                is_member=False,
                explanation="No tetrahydrofolate one-carbon donor/product pair detected",
            )

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.ONE_CARBON_FOLATE_DONORS
            | self.TETRAHYDROFOLATE_FORMS
            | {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.ONE_CARBON_FOLATE_DONORS
            | self.TETRAHYDROFOLATE_FORMS
            | {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        if not left_core or not right_core:
            return ClassificationResult(
                is_member=False,
                explanation="No non-folate donor/acceptor pair to receive the transferred one-carbon unit",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Folate-dependent one-carbon transfer between a tetrahydrofolate carrier and a non-folate acceptor",
        )
