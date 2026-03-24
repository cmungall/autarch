"""Oxidoreductase acting on metal ions.

EC 1.16: Oxidoreductases that alter the oxidation state of metal ions.

These enzymes catalyze metal ion oxidation/reduction:
- Ferroxidase (ceruloplasmin): 4 Fe2+ + O2 + 4H+ → 4 Fe3+ + 2H2O
- Ferric-chelate reductase: Fe3+-chelate + NAD(P)H → Fe2+-chelate + NAD(P)+
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Metal ion ChEBI IDs (oxidized/reduced pairs)
CHEBI_FE2 = "CHEBI:29033"   # iron(2+) / Fe2+ (ferrous)
CHEBI_FE3 = "CHEBI:29034"   # iron(3+) / Fe3+ (ferric)
CHEBI_CU1 = "CHEBI:49552"   # copper(1+) / Cu+
CHEBI_CU2 = "CHEBI:29036"   # copper(2+) / Cu2+

# Metal ion pairs: each tuple is (reduced, oxidized) form
METAL_ION_PAIRS = [
    (CHEBI_FE2, CHEBI_FE3),  # Fe2+/Fe3+
    (CHEBI_CU1, CHEBI_CU2),  # Cu+/Cu2+
]

# All metal ions involved in redox
METAL_IONS = {CHEBI_FE2, CHEBI_FE3, CHEBI_CU1, CHEBI_CU2}


class OxidoreductaseActingOnMetalIons(Oxidoreductase):
    """Catalysis of an oxidation-reduction in which the oxidation state of metal ion is altered.

    EC 1.16.x.x includes:
    - Ferroxidase (ceruloplasmin): 4 Fe2+ + O2 + 4 H+ = 4 Fe3+ + 2 H2O
    - Ferric-chelate reductase: 2 Fe3+-siderophore + NAD(P)H = 2 Fe2+ + NAD(P)+ + siderophore
    - Diferric-transferrin reductase

    Examples:
    - 4 Fe2+ + O2 + 4 H+ = 4 Fe3+ + 2 H2O
    - Fe3+-siderophore + NADPH = Fe2+ + NADP+ + siderophore
    """

    GO_ID = "GO:0016722"  # oxidoreductase activity, acting on metal ions
    EC_NUMBER_PREFIX = "1.16.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on metal ions.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Must have BOTH reduced and oxidized forms of a metal ion on opposite sides
        3. The metal ion must be the PRIMARY substrate, not an electron carrier.
           In EC 1.16, the metal ion IS the donor -- so the only other participants
           should be the electron acceptor (O2, NAD(P)+) and small molecules (H2O, H+).
           If other organic substrates are present (lactate, formate, etc.), the metal
           is acting as an electron carrier for a different EC class.
        """
        # First check if it's an oxidoreductase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        left_chebis = {p.chebi_id for p in reaction.left_participants}
        right_chebis = {p.chebi_id for p in reaction.right_participants}
        all_chebis = left_chebis | right_chebis

        # Check for metal ion redox pairs (both forms must be present on opposite sides)
        has_metal_redox = False
        for reduced, oxidized in METAL_ION_PAIRS:
            if reduced in all_chebis and oxidized in all_chebis:
                if (reduced in left_chebis and oxidized in right_chebis) or \
                   (oxidized in left_chebis and reduced in right_chebis):
                    has_metal_redox = True
                    break

        if not has_metal_redox:
            # Check label as fallback for cobalt, mercury, etc.
            label_lower = reaction.label.lower() if reaction.label else ""
            metal_keywords = ["ferroxidase", "ceruloplasmin", "diferric-transferrin"]
            matched = [kw for kw in metal_keywords if kw in label_lower]
            if matched:
                return ClassificationResult(
                    is_member=True,
                    explanation=f"Metal ion oxidoreductase: label contains {', '.join(matched)}"
                )
            return ClassificationResult(
                is_member=False,
                explanation="No metal ion redox pair (Fe2+/Fe3+, Cu+/Cu2+) detected on opposite sides"
            )

        # Metal redox pair found -- now check if metal is the PRIMARY substrate
        # In EC 1.16, the non-metal participants should only be electron acceptors
        # (O2, NAD(P)+/NAD(P)H) and small molecules (H2O, H+, siderophores).
        # If other organic/inorganic substrates are being oxidized/reduced
        # (lactate, formate, nitrate, sulfite, etc.), the metal is just a carrier.
        SMALL_AND_COFACTOR_CHEBIS = {
            CHEBI_FE2, CHEBI_FE3, CHEBI_CU1, CHEBI_CU2,
            "CHEBI:15377",  # H2O
            "CHEBI:15378",  # H+
            "CHEBI:57540",  # NAD+
            "CHEBI:57945",  # NADH
            "CHEBI:58349",  # NADP+
            "CHEBI:57783",  # NADPH
            "CHEBI:15379",  # O2
            "CHEBI:16240",  # H2O2
            "CHEBI:57692",  # FAD
            "CHEBI:57618",  # FADH2
            "CHEBI:58210",  # FMN
        }

        # Count participants that are NOT metal ions, cofactors, or small molecules
        non_trivial = [
            p for p in reaction.left_participants + reaction.right_participants
            if p.chebi_id not in SMALL_AND_COFACTOR_CHEBIS
        ]

        if len(non_trivial) == 0:
            return ClassificationResult(
                is_member=True,
                explanation="Metal ion oxidoreductase: metal ion is primary substrate (no other organic substrates)"
            )

        # If there are non-trivial participants, check if they look like
        # siderophores/chelates (which ARE part of EC 1.16 ferric-chelate reductase)
        label_lower = reaction.label.lower() if reaction.label else ""
        chelate_terms = ["siderophore", "chelate", "transferrin", "ferritin"]
        if any(term in label_lower for term in chelate_terms):
            return ClassificationResult(
                is_member=True,
                explanation="Metal ion oxidoreductase: ferric-chelate reductase pattern"
            )

        return ClassificationResult(
            is_member=False,
            explanation="Metal ion redox pair present but other substrates indicate metal is electron carrier, not primary substrate"
        )
