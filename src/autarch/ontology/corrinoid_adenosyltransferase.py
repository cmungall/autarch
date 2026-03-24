"""corrinoid adenosyltransferase activity.

Catalysis of the reaction: 2 ATP + 2 corrinoid + reduced [electron-transfer flavoprotein] = 2 adenosylcorrinoid + 3 H+ + oxidized [electron-transfer flavoprotein] + 2 triphosphate. The corrinoid can be cob(II)yrinate a,c diamide, cob(II)inamide or cob(II)alamin.
"""

from collections import Counter

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_ATP, CHEBI_FAD, CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_REDUCED_ETF = "CHEBI:58307"
CHEBI_TRIPHOSPHATE = "CHEBI:18036"

CORRINOID_PRODUCT_PAIRS = {
    "CHEBI:58537": "CHEBI:58503",
    "CHEBI:141013": "CHEBI:2480",
    "CHEBI:16304": "CHEBI:18408",
}


class CorrinoidAdenosyltransferase(ReactionClass):
    """corrinoid adenosyltransferase activity.

    Catalysis of the reaction: 2 ATP + 2 corrinoid + reduced [electron-transfer flavoprotein] = 2 adenosylcorrinoid + 3 H+ + oxidized [electron-transfer flavoprotein] + 2 triphosphate. The corrinoid can be cob(II)yrinate a,c diamide, cob(II)inamide or cob(II)alamin.
    """

    GO_ID = "GO:0008817"
    EC_NUMBER_PREFIX = "2.5.1.17"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_counts = Counter(
            participant.chebi_id for participant in reaction.left_participants if participant.chebi_id
        )
        right_counts = Counter(
            participant.chebi_id for participant in reaction.right_participants if participant.chebi_id
        )

        if left_counts[CHEBI_REDUCED_ETF] != 1 or right_counts[CHEBI_FAD] != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Requires reduced electron-transfer flavoprotein substrate and oxidized flavoprotein product",
            )
        if left_counts[CHEBI_ATP] == 0 or left_counts[CHEBI_ATP] != right_counts[CHEBI_TRIPHOSPHATE]:
            return ClassificationResult(
                is_member=False,
                explanation="ATP consumption must be matched by triphosphate release",
            )

        matched_substrate = next(
            (
                substrate
                for substrate, product in CORRINOID_PRODUCT_PAIRS.items()
                if left_counts[substrate] > 0 and left_counts[substrate] == right_counts[product]
            ),
            None,
        )
        if matched_substrate is None:
            return ClassificationResult(
                is_member=False,
                explanation="No matched corrinoid to adenosylcorrinoid substrate-product pair detected",
            )

        substantive_left = {
            chebi
            for chebi in left_counts
            if chebi not in {matched_substrate, CHEBI_REDUCED_ETF, CHEBI_ATP, CHEBI_H_PLUS}
        }
        substantive_right = {
            chebi
            for chebi in right_counts
            if chebi
            not in {
                CORRINOID_PRODUCT_PAIRS[matched_substrate],
                CHEBI_TRIPHOSPHATE,
                CHEBI_FAD,
                CHEBI_H_PLUS,
            }
        }
        if substantive_left or substantive_right:
            return ClassificationResult(
                is_member=False,
                explanation="Contains additional substantive participants outside corrinoid adenosyltransferase chemistry",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Corrinoid adenosyltransferase: ATP-dependent transfer of adenosyl to a corrinoid with triphosphate release and electron-transfer flavoprotein oxidation",
        )
