"""spermine oxidase (propane-1,3-diamine-forming) activity.

Catalysis of the reaction: H2O + O2 + spermidine/spermine/N(1)-acetylspermine = aldehyde + propane-1,3-diamine + H2O2.
"""

from abc import abstractmethod

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H2O2, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

PRODUCT_PAIRS = {
    "CHEBI:57834": "CHEBI:58264",
    "CHEBI:45725": "CHEBI:58869",
    "CHEBI:58101": "CHEBI:58858",
}
CHEBI_PROPANE_13_DIAMINE = "CHEBI:57484"


class SpermineOxidasePropane13DiamineForming(ReactionClass):
    """Legacy implementation for polyamine oxidase (propane-1,3-diamine-forming).

    Catalysis of the reaction: H2O + O2 + spermidine/spermine/N(1)-acetylspermine = aldehyde + propane-1,3-diamine + H2O2.
    """

    EC_NUMBER_PREFIX = "1.5.3.14"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    @abstractmethod
    def _implementation_only(self) -> None:
        """Mark this legacy implementation class as abstract."""

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebis = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        if CHEBI_O2 not in left_chebis or CHEBI_H2O not in left_chebis:
            return ClassificationResult(is_member=False, explanation="Requires dioxygen and water substrates")
        if CHEBI_H2O2 not in right_chebis or CHEBI_PROPANE_13_DIAMINE not in right_chebis:
            return ClassificationResult(is_member=False, explanation="Requires hydrogen peroxide and propane-1,3-diamine products")

        substrate = next((chebi for chebi in PRODUCT_PAIRS if chebi in left_chebis), None)
        if substrate is None:
            return ClassificationResult(is_member=False, explanation="No supported spermine-oxidase substrate detected")

        aldehyde = PRODUCT_PAIRS[substrate]
        if aldehyde not in right_chebis:
            return ClassificationResult(is_member=False, explanation="Missing the expected aldehyde co-product")

        substantive_left = left_chebis - {CHEBI_O2, CHEBI_H2O, substrate}
        substantive_right = right_chebis - {CHEBI_H2O2, CHEBI_PROPANE_13_DIAMINE, aldehyde}
        if substantive_left or substantive_right:
            return ClassificationResult(is_member=False, explanation="Contains additional substantive participants outside the propane-1,3-diamine-forming branch")

        return ClassificationResult(
            is_member=True,
            explanation="Spermine oxidase (propane-1,3-diamine-forming): oxygen-dependent oxidative cleavage yielding propane-1,3-diamine and an aldehyde co-product",
        )
