"""Dioxygenase reaction classification.

EC 1.13: Oxidoreductases acting on single donors with incorporation of O2.
Dioxygenases incorporate both atoms of molecular oxygen into substrates.

Patterns:

- Substrate + O2 → Product (with both O atoms incorporated)
- Ring-cleaving dioxygenases: catechol + O2 → muconate
- Lipoxygenases: fatty acid + O2 → hydroperoxide
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_O2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_CO2,
)


class OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygen(Oxidoreductase):
    """oxidoreductase acting on single donors with incorporation of molecular oxygen

    EC 1.13.x.x includes:
    - Catechol dioxygenases (ring cleavage)
    - Lipoxygenases (fatty acid oxidation)
    - Alpha-ketoglutarate-dependent dioxygenases
    - Aromatic ring dioxygenases

    Examples:
    - Catechol + O2 → muconate
    - Arachidonate + O2 → hydroperoxyeicosatetraenoate
    """

    GO_ID = "GO:0016701"  # oxidoreductase, acting on single donors with O2 incorporation
    EC_NUMBER_PREFIX = "1.13.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a dioxygenase.

        Strategy:
        1. Must use O2 as substrate
        2. Must NOT produce H2O2 (that would be oxidase, EC 1.1.3)
        3. Look for characteristic ring cleavage or hydroperoxide products
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must use O2 as substrate (ChEBI ID primary)
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        ) or any(x in substrate_str for x in ["dioxygen", "oxygen", " o2"])

        if not has_o2:
            return ClassificationResult(
                is_member=False,
                explanation="No O2 substrate"
            )

        # Dioxygenases do NOT produce H2O2 (oxidases do)
        produces_h2o2 = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        ) or ("peroxide" in product_str and "hydroperoxy" not in product_str)

        if produces_h2o2:
            return ClassificationResult(
                is_member=False,
                explanation="Produces H2O2 - oxidase, not dioxygenase"
            )

        # Monooxygenases produce H2O (one O into substrate, one to H2O)
        # But some dioxygenases also produce H2O with CO2
        produces_h2o = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )

        produces_co2 = any(
            p.chebi_id == CHEBI_CO2
            for p in reaction.right_participants
        )

        # Check for characteristic dioxygenase products

        # Ring cleavage products (muconate, semialdehyde)
        ring_cleavage_patterns = [
            "muconate", "muconic", "semialdehyde",
            "2-hydroxymuconate", "carboxymuconate",
        ]

        has_ring_cleavage = any(rcp in product_str for rcp in ring_cleavage_patterns)

        # Hydroperoxide products (lipoxygenase)
        has_hydroperoxide = "hydroperoxy" in product_str or "hydroperoxide" in product_str

        # Characteristic dioxygenase substrates
        dioxygenase_substrates = [
            "catechol", "protocatechuate", "gentisate", "homogentisate",
            "hydroxybenzoate", "dihydroxybenzoate",
            "tryptophan", "indole",
            "arachidonate", "linolenate", "linoleate", "eicosatetraenoate",
        ]

        has_dioxygenase_substrate = any(ds in substrate_str for ds in dioxygenase_substrates)

        # Ring cleavage dioxygenase
        if has_ring_cleavage:
            return ClassificationResult(
                is_member=True,
                explanation="Dioxygenase: aromatic ring cleavage"
            )

        # Lipoxygenase (fatty acid + O2 → hydroperoxide)
        if has_hydroperoxide:
            return ClassificationResult(
                is_member=True,
                explanation="Dioxygenase: lipoxygenase (hydroperoxide formation)"
            )

        # Characteristic substrate with O2
        if has_dioxygenase_substrate and not produces_h2o2:
            # Exclude if produces H2O without CO2 (likely monooxygenase)
            if produces_h2o and not produces_co2:
                return ClassificationResult(
                    is_member=False,
                    explanation="Produces H2O only - likely monooxygenase"
                )

            return ClassificationResult(
                is_member=True,
                explanation="Dioxygenase: O2 incorporation"
            )

        # Oxidative decarboxylation with O2 (some dioxygenases release CO2)
        if produces_co2 and has_o2:
            # Check for alpha-ketoglutarate dependent pattern
            has_oxoglutarate = "oxoglutarate" in label_lower

            if has_oxoglutarate:
                return ClassificationResult(
                    is_member=False,
                    explanation="Alpha-ketoglutarate dependent - likely EC 1.14"
                )

            # Simple O2 + substrate → products + CO2 pattern
            if not produces_h2o:
                return ClassificationResult(
                    is_member=True,
                    explanation="Dioxygenase: oxidative decarboxylation"
                )

        return ClassificationResult(
            is_member=False,
            explanation="No dioxygenase pattern"
        )
