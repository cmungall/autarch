"""Squalene monooxygenase reaction classifier."""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_O2
from autarch.ontology.oxidoreductase import Oxidoreductase

# Reaction participants for RHEA:25282 / GO:0004506
CHEBI_SQUALENE = "CHEBI:15440"
CHEBI_EPOXYSQUALENE = "CHEBI:15441"
CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE = "CHEBI:57618"
CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE = "CHEBI:58210"


class SqualeneMonooxygenase(Oxidoreductase):
    """squalene monooxygenase"""

    GO_ID = "GO:0004506"  # squalene monooxygenase activity
    EC_NUMBER_PREFIX = "1.14.14.17"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for squalene -> 2,3-epoxysqualene monooxygenation."""
        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        forward = (
            CHEBI_SQUALENE in left_ids
            and CHEBI_EPOXYSQUALENE in right_ids
            and CHEBI_O2 in left_ids
            and CHEBI_H2O in right_ids
        )
        reverse = (
            CHEBI_SQUALENE in right_ids
            and CHEBI_EPOXYSQUALENE in left_ids
            and CHEBI_O2 in right_ids
            and CHEBI_H2O in left_ids
        )

        if not (forward or reverse):
            return ClassificationResult(
                is_member=False,
                explanation="Missing squalene/epoxysqualene monooxygenase signature",
            )

        has_hemoprotein_reductase_cycle = (
            CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE in left_ids
            and CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE in right_ids
        ) or (
            CHEBI_REDUCED_HEMOPROTEIN_REDUCTASE in right_ids
            and CHEBI_OXIDIZED_HEMOPROTEIN_REDUCTASE in left_ids
        )

        if not has_hemoprotein_reductase_cycle:
            return ClassificationResult(
                is_member=False,
                explanation="Missing NADPH--hemoprotein reductase redox pair",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Squalene monooxygenase: squalene epoxidation with O2",
        )
