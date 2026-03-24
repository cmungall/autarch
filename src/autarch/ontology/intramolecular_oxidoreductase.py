"""Intramolecular oxidoreductase reaction classification.

EC 5.3: Intramolecular oxidoreductases (isomerases that catalyze internal redox).
These enzymes catalyze hydrogen transfer within a single molecule.

Patterns:
- Aldose ⟷ Ketose (e.g., glucose-6-phosphate ⟷ fructose-6-phosphate)
- Keto-enol tautomerization
- Double bond migration
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase


class IntramolecularOxidoreductase(Isomerase):
    """intramolecular oxidoreductase

    EC 5.3.x.x includes:
    - Aldose-ketose isomerases (5.3.1): glucose-6-P ⟷ fructose-6-P
    - Keto-enol tautomerases (5.3.2): oxaloacetate (keto) ⟷ enol-oxaloacetate
    - Double bond isomerases (5.3.3): Δ5 → Δ4 steroid isomerase
    - S-S bond isomerases (5.3.4): protein disulfide isomerase

    Examples:
    - D-glucose 6-phosphate = D-fructose 6-phosphate
    - Prostaglandin A1 = prostaglandin C1
    """

    GO_ID = "GO:0016860"  # intramolecular oxidoreductase activity
    EC_NUMBER_PREFIX = "5.3.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an intramolecular oxidoreductase.

        Strategy:
        1. Must be 1 → 1 stoichiometry (isomerization)
        2. No external cofactors (internal H transfer)
        3. Look for characteristic aldose-ketose or keto-enol patterns
        """
        # Must be 1 → 1 (isomerization)
        if len(reaction.left_participants) != 1 or len(reaction.right_participants) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Not 1→1 stoichiometry - not isomerization"
            )

        substrate = reaction.left_participants[0]
        product = reaction.right_participants[0]

        # Same ChEBI ID means no reaction
        if substrate.chebi_id and substrate.chebi_id == product.chebi_id:
            return ClassificationResult(
                is_member=False,
                explanation="Same compound on both sides"
            )

        # Aldose-ketose isomerase patterns (EC 5.3.1)
        aldose_patterns = ["glucose", "mannose", "galactose", "ribose", "arabinose",
                         "xylose", "aldose", "aldohexose", "aldopentose"]
        ketose_patterns = ["fructose", "ribulose", "xylulose", "ketose",
                         "ketohexose", "ketopentose"]

        substrate_name = substrate.name.lower() if substrate.name else ""
        product_name = product.name.lower() if product.name else ""

        has_aldose = any(ap in substrate_name or ap in product_name for ap in aldose_patterns)
        has_ketose = any(kp in substrate_name or kp in product_name for kp in ketose_patterns)

        if has_aldose and has_ketose:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: aldose-ketose isomerase"
            )

        # Sugar-phosphate isomerase (common form of 5.3.1)
        sugar_phosphates = ["phosphate", "6-phosphate", "1-phosphate"]
        (
            any(sp in substrate_name for sp in sugar_phosphates) and
            any(sp in product_name for sp in sugar_phosphates)
        )

        sugar_pair_patterns = [
            ("glucose", "fructose"),
            ("mannose", "fructose"),
            ("ribose", "ribulose"),
            ("xylose", "xylulose"),
            ("arabinose", "ribulose"),
        ]

        is_sugar_isomerization = any(
            (s1 in substrate_name and s2 in product_name) or
            (s2 in substrate_name and s1 in product_name)
            for s1, s2 in sugar_pair_patterns
        )

        if is_sugar_isomerization:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: sugar isomerase"
            )

        # Prostaglandin/eicosanoid double bond isomerization (EC 5.3.3)
        eicosanoid_patterns = ["prostaglandin", "leukotriene", "thromboxane"]
        has_eicosanoid = any(ep in substrate_name or ep in product_name for ep in eicosanoid_patterns)

        if has_eicosanoid:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: eicosanoid isomerase"
            )

        # Carotenoid/terpenoid isomerization (violaxanthin → neoxanthin)
        carotenoid_patterns = ["xanthin", "carotene", "retinal", "lycopene"]
        has_carotenoid = any(cp in substrate_name or cp in product_name for cp in carotenoid_patterns)

        if has_carotenoid:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: carotenoid isomerase"
            )

        # Steroid isomerization (delta-isomerase)
        steroid_patterns = ["steroid", "pregn", "androst", "estro", "cholest"]
        has_steroid = any(sp in substrate_name or sp in product_name for sp in steroid_patterns)

        if has_steroid:
            return ClassificationResult(
                is_member=True,
                explanation="Intramolecular oxidoreductase: steroid isomerase"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No intramolecular oxidoreductase pattern"
        )
