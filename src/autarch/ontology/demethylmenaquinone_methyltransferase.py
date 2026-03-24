"""demethylmenaquinone methyltransferase activity.

Catalysis of the reaction: a 2-demethylmenaquinol + S-adenosyl-L-methionine =
a menaquinol + H(+) + S-adenosyl-L-homocysteine. Reaction substrates can have
varying polyprenyl side chain length.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_SAH, CHEBI_SAM
from autarch.ontology.methyltransferase import Methyltransferase

CHEBI_DEMETHYLMENAQUINONE = "CHEBI:28192"
CHEBI_MENAQUINONE = "CHEBI:16374"
CHEBI_DEMETHYLMENAQUINOL = "CHEBI:55437"
CHEBI_MENAQUINOL = "CHEBI:18151"
CHEBI_DEMETHYLMENAQUINOL_8 = "CHEBI:61873"
CHEBI_MENAQUINOL_8 = "CHEBI:61684"
CHEBI_DEMETHYLMENAQUINOL_7 = "CHEBI:64806"
CHEBI_MENAQUINOL_7 = "CHEBI:64834"
CHEBI_DEMETHYLMENAQUINOL_6 = "CHEBI:84539"
CHEBI_MENAQUINOL_6 = "CHEBI:84536"
CHEBI_DEMETHYLMENAQUINOL_9 = "CHEBI:84542"
CHEBI_MENAQUINOL_9 = "CHEBI:84541"
CHEBI_DEMETHYLMENAQUINOL_10 = "CHEBI:84546"
CHEBI_MENAQUINOL_10 = "CHEBI:84544"
CHEBI_DEMETHYLMENAQUINOL_11 = "CHEBI:84548"
CHEBI_MENAQUINOL_11 = "CHEBI:84547"
CHEBI_DEMETHYLMENAQUINOL_12 = "CHEBI:84551"
CHEBI_MENAQUINOL_12 = "CHEBI:84550"
CHEBI_DEMETHYLMENAQUINOL_13 = "CHEBI:84553"
CHEBI_MENAQUINOL_13 = "CHEBI:84552"
CHEBI_METHOXY_OCTAPRENYLBENZOQUINONE = "CHEBI:28423"
CHEBI_METHYL_METHOXY_OCTAPRENYLBENZOQUINONE = "CHEBI:28636"


class DemethylmenaquinoneMethyltransferase(Methyltransferase):
    """demethylmenaquinone methyltransferase activity.

    Catalysis of the reaction: a 2-demethylmenaquinol +
    S-adenosyl-L-methionine = a menaquinol + H(+) +
    S-adenosyl-L-homocysteine. Reaction substrates can have varying
    polyprenyl side chain length.
    """

    GO_ID = "GO:0043770"
    EC_NUMBER_PREFIX = "2.1.1.163"

    SUBSTRATE_TO_PRODUCT = {
        CHEBI_METHOXY_OCTAPRENYLBENZOQUINONE: CHEBI_METHYL_METHOXY_OCTAPRENYLBENZOQUINONE,
        CHEBI_DEMETHYLMENAQUINONE: CHEBI_MENAQUINONE,
        CHEBI_DEMETHYLMENAQUINOL: CHEBI_MENAQUINOL,
        CHEBI_DEMETHYLMENAQUINOL_6: CHEBI_MENAQUINOL_6,
        CHEBI_DEMETHYLMENAQUINOL_7: CHEBI_MENAQUINOL_7,
        CHEBI_DEMETHYLMENAQUINOL_8: CHEBI_MENAQUINOL_8,
        CHEBI_DEMETHYLMENAQUINOL_9: CHEBI_MENAQUINOL_9,
        CHEBI_DEMETHYLMENAQUINOL_10: CHEBI_MENAQUINOL_10,
        CHEBI_DEMETHYLMENAQUINOL_11: CHEBI_MENAQUINOL_11,
        CHEBI_DEMETHYLMENAQUINOL_12: CHEBI_MENAQUINOL_12,
        CHEBI_DEMETHYLMENAQUINOL_13: CHEBI_MENAQUINOL_13,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return parent_result

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in {CHEBI_SAM, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in {CHEBI_SAH, CHEBI_H_PLUS}
        ]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one demethylmenaquinone substrate and one menaquinone product",
            )

        substrate = left_core[0].chebi_id
        product = right_core[0].chebi_id
        if substrate not in self.SUBSTRATE_TO_PRODUCT:
            return ClassificationResult(
                is_member=False,
                explanation="Requires a demethylmenaquinone or demethylmenaquinol substrate",
            )
        if self.SUBSTRATE_TO_PRODUCT[substrate] != product:
            return ClassificationResult(
                is_member=False,
                explanation="Product does not match the expected methylated menaquinone branch",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Demethylmenaquinone methyltransferase: SAM-dependent methylation of a demethylmenaquinone scaffold",
        )
