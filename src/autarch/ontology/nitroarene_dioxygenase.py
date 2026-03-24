"""nitroarene dioxygenase.

Catalysis of the reaction: a nitroarene + NADH + O2 = a catecholic product +
nitrite + NAD(+).
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.molecules import CHEBI_H_PLUS, CHEBI_NADH, CHEBI_NAD_PLUS, CHEBI_O2
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass

CHEBI_NITRITE = "CHEBI:16301"
NITRO_AROMATIC_PATTERN = Chem.MolFromSmarts("c[N+](=O)[O-]")
CATECHOL_PATTERN = Chem.MolFromSmarts("c([OX2H1,OX1-])c([OX2H1,OX1-])")

KNOWN_NITRO_SUBSTRATES = {
    "CHEBI:27798",
    "CHEBI:957",
    "CHEBI:33098",
    "CHEBI:39931",
    "CHEBI:142281",
    "CHEBI:142283",
    "CHEBI:51398",
    "CHEBI:51397",
    "CHEBI:142285",
    "CHEBI:35227",
    "CHEBI:142287",
    "CHEBI:82420",
}
KNOWN_CATECHOL_PRODUCTS = {
    "CHEBI:18135",
    "CHEBI:142280",
    "CHEBI:18404",
    "CHEBI:17254",
    "CHEBI:142282",
    "CHEBI:142284",
    "CHEBI:57730",
    "CHEBI:142286",
    "CHEBI:142288",
    "CHEBI:27772",
}


class NitroareneDioxygenase(ReactionClass):
    """nitroarene dioxygenase.

    Catalysis of the reaction: a nitroarene + NADH + O2 = a catecholic product
    + nitrite + NAD(+).
    """

    EC_NUMBER_PREFIX = "1.14.12.23"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        left_ids = [participant.chebi_id for participant in reaction.left_participants if participant.chebi_id]
        right_ids = [participant.chebi_id for participant in reaction.right_participants if participant.chebi_id]

        if CHEBI_NADH not in left_ids or CHEBI_O2 not in left_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires NADH and dioxygen cosubstrates",
            )
        if CHEBI_NAD_PLUS not in right_ids or CHEBI_NITRITE not in right_ids:
            return ClassificationResult(
                is_member=False,
                explanation="Requires NAD(+) regeneration and nitrite release",
            )

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in {CHEBI_NADH, CHEBI_O2, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in {CHEBI_NAD_PLUS, CHEBI_NITRITE, CHEBI_H_PLUS}
        ]
        if len(left_core) != 1 or len(right_core) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Expected one nitroarene substrate and one catecholic product",
            )

        substrate = left_core[0]
        product = right_core[0]
        if not self._is_nitroarene(substrate):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate is not a supported nitroarene donor",
            )
        if not self._is_catecholic_product(product):
            return ClassificationResult(
                is_member=False,
                explanation="Product is not a supported catecholic nitroarene dioxygenase product",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Nitroarene dioxygenase: NADH-dependent oxygenation of a nitroarene with nitrite release",
        )

    @staticmethod
    def _is_nitroarene(participant: Participant) -> bool:
        if participant.chebi_id in KNOWN_NITRO_SUBSTRATES:
            return True
        mol = participant.get_mol()
        return mol is not None and NITRO_AROMATIC_PATTERN is not None and mol.HasSubstructMatch(NITRO_AROMATIC_PATTERN)

    @staticmethod
    def _is_catecholic_product(participant: Participant) -> bool:
        if participant.chebi_id in KNOWN_CATECHOL_PRODUCTS:
            return True
        mol = participant.get_mol()
        return mol is not None and CATECHOL_PATTERN is not None and mol.HasSubstructMatch(CATECHOL_PATTERN)
