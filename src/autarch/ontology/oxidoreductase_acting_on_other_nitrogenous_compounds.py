"""Oxidoreductase acting on other nitrogenous compounds as donors.

EC 1.7: Oxidoreductases acting on nitrogenous compounds as donors,
excluding CH-NH2 (EC 1.4) and CH-NH (EC 1.5) groups.

These enzymes catalyze redox reactions on nitrogen oxide species:
- Nitrite oxidase: NO2- → NO3-
- Nitrate reductase: NO3- → NO2-
- Hydroxylamine reductase: NH2OH → NH3
- Nitric oxide reductase: NO → N2O
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Nitrogen oxide species (EC 1.7 donors/products)
CHEBI_NITRITE = "CHEBI:16301"        # NO2-
CHEBI_NITRATE = "CHEBI:17632"        # NO3-
CHEBI_HYDROXYLAMINE = "CHEBI:15429"  # NH2OH
CHEBI_NITRIC_OXIDE = "CHEBI:16480"   # NO
CHEBI_NITROUS_OXIDE = "CHEBI:17045"  # N2O
CHEBI_DINITROGEN = "CHEBI:17997"     # N2

NITROGEN_OXIDE_SPECIES = {
    CHEBI_NITRITE,
    CHEBI_NITRATE,
    CHEBI_HYDROXYLAMINE,
    CHEBI_NITRIC_OXIDE,
    CHEBI_NITROUS_OXIDE,
    CHEBI_DINITROGEN,
}


class OxidoreductaseActingOnOtherNitrogenousCompounds(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which a nitrogenous group, excluding NH and NH2 groups, acts as a hydrogen or electron donor and reduces a hydrogen or electron acceptor.

    EC 1.7.x.x includes:
    - Nitrite oxidase: nitrite + acceptor = nitrate + reduced acceptor
    - Nitrate reductase: nitrate + donor = nitrite + oxidized donor
    - Hydroxylamine reductase: hydroxylamine + acceptor = NH3 + oxidized acceptor
    - Nitric oxide reductase: NO + donor = N2O + oxidized donor

    Examples:
    - nitrite + NAD+ + H2O = nitrate + NADH + H+
    - hydroxylamine + acceptor = NH3 + H2O + oxidized acceptor
    - 2 NO + 2 ferrocytochrome c = N2O + H2O + 2 ferricytochrome c
    """

    GO_ID = "GO:0016661"  # oxidoreductase activity, acting on other nitrogenous compounds as donors
    EC_NUMBER_PREFIX = "1.7.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on nitrogenous compounds.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Must involve nitrogen oxide species (nitrite, nitrate, hydroxylamine, NO, N2O)
        3. Exclude reactions where ammonia/ammonium is the only nitrogen species
           (those are CH-NH2 oxidoreductases, EC 1.4)
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        all_participants = reaction.left_participants + reaction.right_participants

        # Check for nitrogen oxide species
        nitrogen_species_found = [
            p.chebi_id for p in all_participants
            if p.chebi_id in NITROGEN_OXIDE_SPECIES
        ]

        if nitrogen_species_found:
            return ClassificationResult(
                is_member=True,
                explanation=f"Nitrogenous compound oxidoreductase: nitrogen oxide species present ({', '.join(nitrogen_species_found)})"
            )

        # Also check reaction label for nitrogen compound keywords
        label_lower = reaction.label.lower() if reaction.label else ""
        nitrogen_keywords = ["nitrite", "nitrate", "hydroxylamine", "nitric oxide",
                            "nitrous oxide", "nitrogenase", "nitrogen fixation"]
        matched_keywords = [kw for kw in nitrogen_keywords if kw in label_lower]

        if matched_keywords:
            return ClassificationResult(
                is_member=True,
                explanation=f"Nitrogenous compound oxidoreductase: label contains {', '.join(matched_keywords)}"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No nitrogen oxide species (nitrite, nitrate, hydroxylamine, NO, N2O) detected"
        )
