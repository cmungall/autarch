"""tRNA-dihydrouridine(16/17) synthase [NAD(P)(+)].

Catalysis of reversible NAD(P)-linked interconversion between tRNA
dihydrouridine and uridine states for the 16/17 branch.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.gdp_4_dehydro_d_rhamnose_reductase import _check_reversible_redox_pair
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_TRNA_DIHYDROURIDINE = "CHEBI:74443"
CHEBI_TRNA_URIDINE = "CHEBI:65315"


class TRNADihydrouridine1617SynthaseNADP(ReactionClass):
    """tRNA-dihydrouridine(16/17) synthase [NAD(P)(+)].

    Catalysis of reversible NAD(P)-linked interconversion between tRNA
    dihydrouridine and uridine states for the 16/17 branch.

    The current structured RHEA cache represents the modified tRNA residues
    generically rather than position-specifically, so this classifier is an
    intentionally coarse structural approximation of the 16/17 branch.
    """

    EC_NUMBER_PREFIX = "1.3.1.88"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return _check_reversible_trna_dihydrouridine_pair(
            reaction,
            explanation="tRNA-dihydrouridine(16/17) synthase [NAD(P)(+)]: reversible NAD(P)-linked tRNA dihydrouridine interconversion",
        )


def _check_reversible_trna_dihydrouridine_pair(
    reaction: Reaction,
    explanation: str,
) -> ClassificationResult:
    return _check_reversible_redox_pair(
        reaction,
        substrate_to_product={CHEBI_TRNA_DIHYDROURIDINE: CHEBI_TRNA_URIDINE},
        cofactor_pairs={CHEBI_NADP_PLUS: CHEBI_NADPH, CHEBI_NAD_PLUS: CHEBI_NADH},
        explanation=explanation,
    )
