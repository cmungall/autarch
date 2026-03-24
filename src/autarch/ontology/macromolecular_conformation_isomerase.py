"""Macromolecular conformation isomerase reaction classification.

EC 5.6: Isomerases that alter macromolecular conformation.
These enzymes change the topology or conformation of DNA, RNA,
or proteins without changing covalent bonds (or by transiently
breaking and rejoining them).
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase


class MacromolecularConformationIsomerase(Isomerase):
    """isomerase altering macromolecular conformation

    EC 5.6.x.x includes:
    - DNA topoisomerases (5.6.2): change DNA supercoiling
    - RNA topoisomerases
    - Chaperones/chaperonins that alter protein conformation

    Examples:
    - DNA (supercoiled) = DNA (relaxed) via topoisomerase
    - Misfolded protein = correctly folded protein via chaperone
    """

    GO_ID = "GO:0120543"  # isomerase activity, altering macromolecular conformation
    EC_NUMBER_PREFIX = "5.6.-.-"

    CONFORMATION_LABEL_PATTERNS = [
        "topoisomerase", "gyrase",
        "chaperone", "chaperonin",
        "helicase",  # some helicases are EC 5.6
        "supercoil", "relaxase",
        "recombinase",
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a macromolecular conformation isomerase.

        Strategy:
        1. These reactions involve macromolecular conformation changes
        2. Label-based detection since the substrates are macromolecules
           not easily represented as SMILES
        """
        label_lower = reaction.label.lower() if reaction.label else ""

        has_conformation_pattern = any(
            pattern in label_lower
            for pattern in self.CONFORMATION_LABEL_PATTERNS
        )

        if has_conformation_pattern:
            return ClassificationResult(
                is_member=True,
                explanation="Macromolecular conformation isomerase: topology/conformation change"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No macromolecular conformation isomerase pattern"
        )
