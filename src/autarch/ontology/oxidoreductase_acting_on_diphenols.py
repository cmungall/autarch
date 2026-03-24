"""Oxidoreductase acting on diphenols and related substances as donors.

EC 1.10: Oxidoreductases acting on diphenols and related substances.

These enzymes use quinol/diphenol substrates as electron donors:
- Laccase: diphenol + O2 → quinone + H2O
- Catechol oxidase: catechol + O2 → o-quinone + H2O
- Quinol oxidase: ubiquinol + O2 → ubiquinone + H2O
- Cytochrome bc1 complex: ubiquinol + cytochrome c
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Quinol/quinone pairs (EC 1.10 substrates)
CHEBI_UBIQUINOL = "CHEBI:17976"       # reduced ubiquinone (QH2)
CHEBI_UBIQUINONE = "CHEBI:16389"      # oxidized ubiquinone (Q)
CHEBI_PLASTOQUINOL = "CHEBI:28108"    # reduced plastoquinone (PQH2)
CHEBI_PLASTOQUINONE = "CHEBI:27594"   # oxidized plastoquinone (PQ)
CHEBI_MENAQUINOL = "CHEBI:18151"      # reduced menaquinone
CHEBI_MENAQUINONE = "CHEBI:18151"     # menaquinone (vitamin K2)
CHEBI_CATECHOL = "CHEBI:18135"        # catechol (1,2-benzenediol)
CHEBI_ASCORBATE = "CHEBI:29073"       # ascorbate (vitamin C)
CHEBI_DEHYDROASCORBATE = "CHEBI:17242"  # dehydroascorbate

# All diphenol/quinol donor species
DIPHENOL_DONORS = {
    CHEBI_UBIQUINOL,
    CHEBI_UBIQUINONE,
    CHEBI_PLASTOQUINOL,
    CHEBI_PLASTOQUINONE,
    CHEBI_CATECHOL,
    CHEBI_ASCORBATE,
    CHEBI_DEHYDROASCORBATE,
}


class OxidoreductaseActingOnDiphenols(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which a diphenol or related substance acts as a hydrogen or electron donor and reduces a hydrogen or electron acceptor.

    EC 1.10.x.x includes:
    - Laccase: 4 benzenediol + O2 = 4 benzosemiquinone + 2 H2O
    - Catechol oxidase: 2 catechol + O2 = 2 o-quinone + 2 H2O
    - Ubiquinol-cytochrome c reductase (Complex III)
    - Plastoquinol-plastocyanin reductase
    - Ascorbate oxidase

    Examples:
    - ubiquinol + 2 ferricytochrome c = ubiquinone + 2 ferrocytochrome c + 2 H+
    - 2 catechol + O2 = 2 1,2-benzoquinone + 2 H2O
    - 2 L-ascorbate + O2 = 2 dehydroascorbate + 2 H2O
    """

    GO_ID = "GO:0016679"  # oxidoreductase activity, acting on diphenols and related substances as donors
    EC_NUMBER_PREFIX = "1.10.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on diphenols.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Must have diphenol/quinol species as DONORS (reactants on left side)
           EC 1.10 specifically requires diphenol as the electron donor
        3. Quinol/quinone pairs: ubiquinol/ubiquinone, plastoquinol/plastoquinone
        4. Catechol as reactant (not product -- product catechol indicates different EC class)
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Quinol/quinone pairs -- presence of EITHER form indicates EC 1.10
        # because direction-neutral RHEA can have either as reactant
        QUINOL_QUINONE_SPECIES = {
            CHEBI_UBIQUINOL, CHEBI_UBIQUINONE,
            CHEBI_PLASTOQUINOL, CHEBI_PLASTOQUINONE,
        }

        all_participants = reaction.left_participants + reaction.right_participants

        # Check for quinol/quinone species (these are definitive for EC 1.10)
        quinol_species = [
            p.chebi_id for p in all_participants
            if p.chebi_id in QUINOL_QUINONE_SPECIES
        ]

        if quinol_species:
            return ClassificationResult(
                is_member=True,
                explanation=f"Diphenol oxidoreductase: quinol/quinone species present ({', '.join(quinol_species)})"
            )

        # Check for catechol as REACTANT (donor), not as product
        # Catechol as product (e.g., salicylate → catechol) is a monooxygenase, not EC 1.10
        has_catechol_reactant = any(
            p.chebi_id == CHEBI_CATECHOL
            for p in reaction.left_participants
        )
        has_catechol_product = any(
            p.chebi_id == CHEBI_CATECHOL
            for p in reaction.right_participants
        )

        if has_catechol_reactant and not has_catechol_product:
            return ClassificationResult(
                is_member=True,
                explanation="Diphenol oxidoreductase: catechol as electron donor (reactant)"
            )

        # Check label for specific quinol/quinone terms (not substring-prone ones)
        label_lower = reaction.label.lower() if reaction.label else ""
        # Use word-boundary-safe keywords: "ubiquinol" not "quinol" (avoids "quinoline")
        diphenol_keywords = [
            "ubiquinol", "ubiquinone", "plastoquinol", "plastoquinone",
            "menaquinol", "menaquinone", "diphenol", "hydroquinone",
            "benzenediol",
        ]
        matched_keywords = [kw for kw in diphenol_keywords if kw in label_lower]

        if matched_keywords:
            return ClassificationResult(
                is_member=True,
                explanation=f"Diphenol oxidoreductase: label contains {', '.join(matched_keywords)}"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No diphenol/quinol donor species detected"
        )
