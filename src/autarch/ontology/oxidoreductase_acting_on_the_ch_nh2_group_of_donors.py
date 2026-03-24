"""Flavin oxidase reaction classification using pattern DSL.

Oxidoreductases acting on CH-NH groups (GO:0016645).
This includes:
- Flavin oxidases (O2 → H2O2): amino acid oxidases, glucose oxidase
- NAD(P)-linked dehydrogenases: proline dehydrogenase, saccharopine dehydrogenase
- Electron-transfer flavoprotein-linked reactions

History
-------

## 2025-12-22

Fixed GO_ID mismatch: was using GO:0016645 (CH-NH, secondary amines, EC 1.5) but
EC_NUMBER_PREFIX was 1.4.3.- (CH-NH2, primary amines). Changed to GO:0016638
which correctly matches EC 1.4 (amino acid oxidases, amine oxidases).

Previously fixed inheritance: was Oxidase (requires O2) but GO:0016645 includes
NAD(P)-linked reactions that don't use O2. Changed to Oxidoreductase base class.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.pattern_dsl import var
from autarch.molecules import (
    CHEBI_H2O2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    p,
    water,
)


class OxidoreductaseActingOnTheCHNH2GroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on the CH-NH2 group of donors

    Examples:
    - D-amino acid oxidase: D-amino acid + O2 + H2O → α-keto acid + NH3 + H2O2
    - L-amino acid oxidase: L-amino acid + O2 + H2O → α-keto acid + NH3 + H2O2
    - Glucose oxidase: glucose + O2 → gluconic acid + H2O2
    - Cholesterol oxidase: cholesterol + O2 → cholest-4-en-3-one + H2O2
    """

    GO_ID = "GO:0016638"  # oxidoreductase activity, acting on the CH-NH2 group of donors
    EC_NUMBER_PREFIX = "1.4.-.-"  # oxidoreductases acting on CH-NH2 group

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Amino acid oxidase: amino acid + O2 + H2O → keto acid + NH3 + H2O2
        var("amino_acid") + var("oxygen") + p(water) 
        >> var("keto_acid") + var("ammonia") + var("h2o2"),
        # General flavin oxidase: substrate + O2 → oxidized substrate + H2O2
        var("substrate") + var("oxygen") >> var("oxidized_substrate") + var("h2o2"),
        # With water involvement: substrate + O2 + H2O → product + H2O2 + byproduct
        var("substrate") + var("oxygen") + p(water) 
        >> var("product") + var("h2o2") + var("byproduct"),
        # Glucose oxidase type: sugar + O2 → acid + H2O2
        var("sugar") + var("oxygen") >> var("acid") + var("h2o2"),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on CH-NH groups.

        GO:0016645 is "oxidoreductase activity, acting on the CH-NH group of donors".
        This is a broad class including:
        1. Flavin oxidases: O2 → H2O2 (amino acid oxidases, glucose oxidase, etc.)
        2. NAD(P)-linked dehydrogenases acting on CH-NH groups

        Strategy:
        - Flavin oxidase path: O2 consumption AND H2O2 production (primary indicator)
        - Dehydrogenase path: NAD(P)+ consumption AND NAD(P)H production AND ammonia release
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        # Check for O2 consumption (flavin oxidase path)
        has_oxygen = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )

        # Check for H2O2 production (flavin oxidase signature)
        has_h2o2_product = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        )

        # Exclude peroxidases (use H2O2 as substrate)
        has_h2o2_substrate = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.left_participants
        )

        if has_h2o2_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Uses H2O2 as substrate - peroxidase, not oxidase"
            )

        # Flavin oxidase path: O2 + substrate → H2O2 + product
        is_flavin_oxidase = has_oxygen and has_h2o2_product

        if is_flavin_oxidase:
            # Check for ammonia production (required for CH-NH2 oxidation)
            # GO:0016638 specifically requires "acting on the CH-NH2 group"
            # Non-amine oxidases (glucose oxidase, alcohol oxidase) don't produce ammonia
            ammonia_ids = {CHEBI_NH3, CHEBI_NH4}
            has_ammonia = any(
                p.chebi_id in ammonia_ids
                for p in reaction.right_participants
            )

            if has_ammonia:
                return ClassificationResult(
                    is_member=True,
                    explanation="CH-NH2 oxidoreductase: O2→H2O2 + ammonia (amino/amine oxidase)"
                )
            # Without ammonia, this is likely a different type of oxidase (glucose, alcohol, etc.)
            # Don't classify as CH-NH2 oxidoreductase

        # Dehydrogenase path: NAD(P)+ as acceptor, produces ammonia
        # (More restricted to avoid FPs - must release ammonia to confirm CH-NH oxidation)
        has_nad_acceptor = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
            for p in reaction.left_participants
        )

        has_nadh_product = any(
            p.chebi_id in {CHEBI_NADH, CHEBI_NADPH}
            for p in reaction.right_participants
        )

        ammonia_ids = {CHEBI_NH3, "CHEBI:28938"}
        has_ammonia = any(
            p.chebi_id in ammonia_ids
            for p in reaction.right_participants
        )

        # For NAD(P)-linked path, require ammonia production to confirm CH-NH oxidation
        if has_nad_acceptor and has_nadh_product and has_ammonia:
            return ClassificationResult(
                is_member=True,
                explanation="CH-NH oxidoreductase: NAD(P)-linked dehydrogenase [produces ammonia]"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No O2→H2O2 or NAD(P)+→NAD(P)H+NH3 pattern"
        )
