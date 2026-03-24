"""L-glutamate gamma-semialdehyde dehydrogenase.

Catalysis of oxidation of L-glutamate gamma-semialdehyde or related cyclic
imines to glutamate derivatives with NAD(P)(+) as acceptor.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.molecules import CHEBI_H2O, CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NADPH, CHEBI_NADP_PLUS, CHEBI_NAD_PLUS
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass


class LGlutamateGammaSemialdehydeDehydrogenase(ReactionClass):
    """L-glutamate gamma-semialdehyde dehydrogenase.

    Catalysis of oxidation of L-glutamate gamma-semialdehyde or related cyclic
    imines to glutamate derivatives with NAD(P)(+) as acceptor.
    """

    EC_NUMBER_PREFIX = "1.2.1.88"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    SUBSTRATE_TO_PRODUCT = {
        "CHEBI:17388": "CHEBI:29985",
        "CHEBI:62612": "CHEBI:6331",
        "CHEBI:58066": "CHEBI:29985",
    }
    COFACTOR_PAIRS = {
        CHEBI_NAD_PLUS: CHEBI_NADH,
        CHEBI_NADP_PLUS: CHEBI_NADPH,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if left_ids.count(CHEBI_H2O) not in {1, 2}:
            return ClassificationResult(is_member=False, explanation="Requires one or two water molecules as reactants")

        matched_pair = next(((left, right) for left, right in self.COFACTOR_PAIRS.items() if left in left_ids and right in right_ids), None)
        if matched_pair is None:
            return ClassificationResult(is_member=False, explanation="Requires NAD(P)(+) reduction to NAD(P)H")

        left_core = [chebi_id for chebi_id in left_ids if chebi_id not in {matched_pair[0], CHEBI_H2O}]
        right_core = [chebi_id for chebi_id in right_ids if chebi_id not in {matched_pair[1], CHEBI_H_PLUS}]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(is_member=False, explanation="Expected one aldehyde/imine substrate and one glutamate-like product")

        substrate = left_core[0]
        product = right_core[0]
        if self.SUBSTRATE_TO_PRODUCT.get(substrate) != product:
            return ClassificationResult(is_member=False, explanation="Organic branch does not match a supported glutamate gamma-semialdehyde dehydrogenase pair")

        return ClassificationResult(is_member=True, explanation="L-glutamate gamma-semialdehyde dehydrogenase: NAD(P)-dependent oxidation of glutamate semialdehyde or a related cyclic imine")
