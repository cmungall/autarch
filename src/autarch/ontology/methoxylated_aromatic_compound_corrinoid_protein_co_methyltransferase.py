"""methoxylated aromatic compound--corrinoid protein Co-methyltransferase.

Catalysis of methyl transfer between a methoxylated aromatic compound and a
corrinoid protein cobalt center.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H_PLUS
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_CORRINOID_COI = "CHEBI:85033"
CHEBI_CORRINOID_METHYL_COIII = "CHEBI:85035"


class MethoxylatedAromaticCompoundCorrinoidProteinCoMethyltransferase(ReactionClass):
    """methoxylated aromatic compound--corrinoid protein Co-methyltransferase.

    Catalysis of methyl transfer between a methoxylated aromatic compound and a
    corrinoid protein cobalt center.
    """

    EC_NUMBER_PREFIX = "2.1.1.382"
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE

    DEMETHYL_TO_METHOXY = {
        "CHEBI:132111": "CHEBI:19950",
        "CHEBI:28591": "CHEBI:18135",
        "CHEBI:63797": "CHEBI:36241",
        "CHEBI:52678": "CHEBI:27810",
        "CHEBI:59128": "CHEBI:30762",
        "CHEBI:192244": "CHEBI:33853",
        "CHEBI:16632": "CHEBI:36241",
    }
    METHOXY_TO_DEMETHYL = {product: substrate for substrate, product in DEMETHYL_TO_METHOXY.items()}
    CORRINOID_TO_AROMATIC_METHYLATION = {
        "CHEBI:52678": "CHEBI:59114",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        forward = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            substrate_to_product=self.DEMETHYL_TO_METHOXY,
            donor=CHEBI_CORRINOID_COI,
            acceptor=CHEBI_CORRINOID_METHYL_COIII,
            donor_requires_hydron=True,
            acceptor_requires_hydron=False,
            explanation="Methoxylated aromatic compound--corrinoid protein Co-methyltransferase: methyl transfer from a methoxylated aromatic compound to corrinoid protein cobalt",
        )
        if forward.is_member:
            return forward

        reverse = self._check_direction(
            left_ids=[participant.chebi_id for participant in reaction.left_participants if participant.chebi_id],
            right_ids=[participant.chebi_id for participant in reaction.right_participants if participant.chebi_id],
            substrate_to_product=self.CORRINOID_TO_AROMATIC_METHYLATION,
            donor=CHEBI_CORRINOID_METHYL_COIII,
            acceptor=CHEBI_CORRINOID_COI,
            donor_requires_hydron=False,
            acceptor_requires_hydron=True,
            explanation="Methoxylated aromatic compound--corrinoid protein Co-methyltransferase: reverse reaction orientation",
        )
        if reverse.is_member:
            return reverse

        return forward

    def _check_direction(
        self,
        left_ids: list[str],
        right_ids: list[str],
        substrate_to_product: dict[str, str],
        donor: str,
        acceptor: str,
        donor_requires_hydron: bool,
        acceptor_requires_hydron: bool,
        explanation: str,
    ) -> ClassificationResult:
        if donor not in left_ids or acceptor not in right_ids:
            return ClassificationResult(is_member=False, explanation="Missing corrinoid methyl-transfer cofactor transition")
        if donor_requires_hydron and CHEBI_H_PLUS not in left_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron with Co(I) corrinoid substrate")
        if acceptor_requires_hydron and CHEBI_H_PLUS not in right_ids:
            return ClassificationResult(is_member=False, explanation="Requires hydron with Co(I) corrinoid product")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {donor, CHEBI_H_PLUS}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {acceptor, CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one aromatic substrate and one aromatic product")

        substrate = left_core[0]
        product = right_core[0]
        if substrate_to_product.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Aromatic branch does not match a supported corrinoid methyl-transfer pair")

        return ClassificationResult(is_member=True, explanation=explanation)
