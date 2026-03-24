"""pyridoxamine-phosphate oxidase activity.

Catalysis of the reaction: pyridoxamine 5'-phosphate + O2 + H2O = pyridoxal 5'-phosphate + H2O2 + NH4(+).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H2O2, CHEBI_NH4, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_PYRIDOXAMINE_PHOSPHATE = "CHEBI:58451"
CHEBI_PYRIDOXAL_PHOSPHATE = "CHEBI:597326"


class PyridoxaminePhosphateOxidase(ReactionClass):
    """pyridoxamine-phosphate oxidase activity.

    Catalysis of the reaction: pyridoxamine 5'-phosphate + O2 + H2O = pyridoxal 5'-phosphate + H2O2 + NH4(+).
    """

    GO_ID = "GO:0004733"
    EC_NUMBER_PREFIX = "1.4.3.5"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_chebis = {participant.chebi_id for participant in reaction.left_participants if participant.chebi_id}
        right_chebis = {participant.chebi_id for participant in reaction.right_participants if participant.chebi_id}

        if left_chebis != {CHEBI_PYRIDOXAMINE_PHOSPHATE, CHEBI_O2, CHEBI_H2O}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires pyridoxamine 5'-phosphate, dioxygen, and water as the complete substrate set",
            )
        if right_chebis != {CHEBI_PYRIDOXAL_PHOSPHATE, CHEBI_H2O2, CHEBI_NH4}:
            return ClassificationResult(
                is_member=False,
                explanation="Requires pyridoxal 5'-phosphate, hydrogen peroxide, and ammonium as the complete product set",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Pyridoxamine phosphate oxidase: oxygen-dependent oxidation of pyridoxamine phosphate to pyridoxal phosphate with ammonium and hydrogen peroxide release",
        )
