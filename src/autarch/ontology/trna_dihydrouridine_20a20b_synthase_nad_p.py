"""tRNA-dihydrouridine(20a/20b) synthase [NAD(P)(+)].

Catalysis of reversible NAD(P)-linked interconversion between tRNA
dihydrouridine and uridine states for the 20a/20b branch.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.ontology.trna_dihydrouridine_1617_synthase_nad_p import (
    _check_reversible_trna_dihydrouridine_pair,
)


class TRNADihydrouridine20A20BSynthaseNADP(ReactionClass):
    """tRNA-dihydrouridine(20a/20b) synthase [NAD(P)(+)].

    Catalysis of reversible NAD(P)-linked interconversion between tRNA
    dihydrouridine and uridine states for the 20a/20b branch.

    The current structured RHEA cache represents the modified tRNA residues
    generically rather than position-specifically, so this classifier is an
    intentionally coarse structural approximation of the 20a/20b branch.
    """

    EC_NUMBER_PREFIX = "1.3.1.90"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_trna_dihydrouridine_pair(
            reaction,
            explanation="tRNA-dihydrouridine(20a/20b) synthase [NAD(P)(+)]: reversible NAD(P)-linked tRNA dihydrouridine interconversion",
        )
