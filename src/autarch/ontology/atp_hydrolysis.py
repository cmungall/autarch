"""ATPase reaction classification using pattern DSL.

ATPases are hydrolases that specifically hydrolyze ATP to provide energy.

NOTE: Current evaluation dataset contains ZERO ATPase reactions (GO:0016887).
This suggests ATPases may not be represented in current RHEA structural data,
or transport ATPases are classified differently than expected.

Previous implementation had 7 false positives with 0 true positives.
Simplified to avoid misclassifying transport reactions.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_PHOSPHATE,
    adp,
    atp,
    h_plus,
    p,
    phosphate,
    water,
)


class ATPHydrolysis(Hydrolase):
    """ATP hydrolysis
    Catalysis of the reaction: ATP + H2O = ADP + H+ + phosphate. ATP hydrolysis is used in some reactions as an energy source, for example to catalyze a reaction or drive transport against a concentration gradient.
    ATP + H2O → AMP + PPi + H+

    Examples:
    - F1F0-ATPase (ATP synthase in reverse)
    - Myosin ATPase (muscle contraction)
    - Na+/K+-ATPase (ion transport)
    - Helicases (DNA unwinding with ATP)
    """

    GO_ID = "GO:0016887"  # ATP hydrolysis activity
    EC_NUMBER_PREFIX = "3.6.-.-"
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE

    # Define patterns using DSL with operator overloading
    PATTERNS: list[Reaction] = [
        # ATP hydrolysis to ADP: ATP + H2O → ADP + Pi + H+?
        p(atp) + p(water) >> p(adp) + p(phosphate) + optional(h_plus),
        # With additional substrates (e.g., coupled reactions)
        # ATP + H2O + substrate → ADP + Pi + product + H+?
        p(atp) + p(water) + var("substrate")
        >> p(adp) + p(phosphate) + var("product") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an ATPase using SIMPLE approach.

        Classify only direct ATP hydrolysis reactions.
        """
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport-coupled ATP usage is classified under transport terms"
            )

        left_ids = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_ids = {p.chebi_id for p in reaction.right_participants if p.chebi_id}

        has_atp = CHEBI_ATP in left_ids
        has_water = CHEBI_H2O in left_ids
        if not (has_atp and has_water):
            return ClassificationResult(
                is_member=False,
                explanation="Missing ATP + H2O substrates"
            )

        allowed_left = {CHEBI_ATP, CHEBI_H2O, CHEBI_H_PLUS}
        if any(cid not in allowed_left for cid in left_ids):
            return ClassificationResult(
                is_member=False,
                explanation="Coupled ATP reaction (extra substrates) - not pure ATP hydrolysis"
            )

        has_phosphate = bool({CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367"} & right_ids)
        has_adp_pi = CHEBI_ADP in right_ids and has_phosphate

        allowed_right = {CHEBI_ADP, CHEBI_PHOSPHATE, "CHEBI:16838", "CHEBI:18367", CHEBI_H_PLUS}
        if not has_adp_pi or any(cid not in allowed_right for cid in right_ids):
            return ClassificationResult(
                is_member=False,
                explanation="Products do not match ATP -> ADP + phosphate hydrolysis signature"
            )

        return ClassificationResult(
            is_member=True,
            explanation="ATP hydrolysis: ATP + H2O -> ADP + phosphate"
        )
