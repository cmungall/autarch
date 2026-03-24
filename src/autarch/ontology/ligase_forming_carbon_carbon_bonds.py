"""Ligase forming carbon-carbon bonds classification.

EC 6.4 ligases catalyze ATP-dependent carboxylation reactions,
incorporating CO2/bicarbonate into organic substrates to form C-C bonds.
Most are biotin-dependent carboxylases.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ligase import Ligase
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_COA,
)

# Bicarbonate (hydrogencarbonate)
CHEBI_BICARBONATE = "CHEBI:17544"


class LigaseFormingCarbonCarbonBonds(Ligase):
    """ligase activity, forming carbon-carbon bonds

    ATP-dependent carboxylation reactions that incorporate CO2 or bicarbonate
    into organic substrates to form new C-C bonds.

    Examples:
    - Pyruvate carboxylase: pyruvate + CO2 + ATP -> oxaloacetate + ADP + Pi
    - Acetyl-CoA carboxylase: acetyl-CoA + CO2 + ATP -> malonyl-CoA + ADP + Pi
    - Propionyl-CoA carboxylase: propionyl-CoA + CO2 + ATP -> methylmalonyl-CoA + ADP + Pi

    Most EC 6.4 enzymes are biotin-dependent carboxylases.
    """

    GO_ID = "GO:0016885"  # ligase activity, forming carbon-carbon bonds
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "6.4.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a C-C bond-forming ligase.

        Strategy:
        1. Must be a ligase (ATP-dependent bond formation)
        2. Must have CO2 or bicarbonate as reactant (carboxylation)
        3. Label-based detection for biotin/carboxylase terms
        """
        # First check if it's a ligase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a ligase: {parent_result.explanation}",
            )

        # Check for CO2 or bicarbonate as reactant
        has_co2 = any(
            p.chebi_id in {CHEBI_CO2, CHEBI_BICARBONATE}
            for p in reaction.left_participants
        )

        if has_co2:
            # Check if CoA substrate is involved (CoA carboxylases)
            has_coa = any(
                p.chebi_id == CHEBI_COA
                or (p.name and "coa" in p.name.lower())
                for p in reaction.left_participants
            )
            substrate_type = "CoA carboxylation" if has_coa else "carboxylation"
            return ClassificationResult(
                is_member=True,
                explanation=f"C-C ligase: ATP-dependent {substrate_type} (CO2 fixation)",
            )

        # Label-based detection for carboxylases and biotin enzymes
        all_names = [
            p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
            if p.name
        ]
        label_terms = {"carboxylase", "biotin"}
        if any(term in name for name in all_names for term in label_terms):
            return ClassificationResult(
                is_member=True,
                explanation="C-C ligase: carboxylase/biotin-dependent enzyme detected by label",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No CO2/bicarbonate substrate - not a C-C bond-forming ligase",
        )
