"""Sulfur oxidoreductase reaction classification.

EC 1.8: Oxidoreductases acting on a sulfur group of donors.
These enzymes catalyze oxidation/reduction of sulfur-containing compounds.

Patterns:
- R-SH + acceptor = R-S-S-R + reduced acceptor (disulfide formation)
- R-S-S-R + donor = 2 R-SH + oxidized donor (disulfide reduction)
- Sulfite + acceptor = Sulfate + reduced acceptor
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)

# Sulfur-containing molecules
CHEBI_GLUTATHIONE = "CHEBI:131728"  # GSH - reduced
CHEBI_GLUTATHIONE_DISULFIDE = "CHEBI:58297"  # GSSG - oxidized
CHEBI_THIOREDOXIN_DITHIOL = "CHEBI:29950"  # reduced thioredoxin
CHEBI_THIOREDOXIN_DISULFIDE = "CHEBI:50058"  # oxidized thioredoxin
CHEBI_CYSTEINE = "CHEBI:35235"
CHEBI_CYSTINE = "CHEBI:17376"  # oxidized cysteine disulfide
CHEBI_SULFITE = "CHEBI:17359"
CHEBI_SULFATE = "CHEBI:16189"
CHEBI_HYDROGEN_SULFIDE = "CHEBI:29919"


class OxidoreductaseActingOnASulfurGroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on a sulfur group of donors

    EC 1.8.x.x includes:
    - Thioredoxin reductase (NADPH + thioredoxin disulfide → thioredoxin + NADP+)
    - Glutathione reductase (NADPH + GSSG → 2 GSH + NADP+)
    - Sulfite reductase (sulfite + NADPH → H2S + NADP+)
    - Disulfide oxidoreductases

    Examples:
    - 2 glutathione + NADP+ = glutathione disulfide + NADPH + H+
    - Thioredoxin + NADP+ = thioredoxin disulfide + NADPH + H+
    """

    GO_ID = "GO:0016667"  # oxidoreductase activity, acting on a sulfur group
    EC_NUMBER_PREFIX = "1.8.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a sulfur oxidoreductase.

        Strategy:
        1. Must be an oxidoreductase (have redox cofactor)
        2. Look for sulfur-containing substrates/products
        3. Pattern: thiol/disulfide interconversion or sulfite/sulfate redox
        """
        # Check for redox cofactors (but exclude iron-sulfur clusters used structurally)
        redox_cofactors = {
            CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
            "CHEBI:30753",  # cytochrome c
        }

        has_redox = any(
            p.chebi_id in redox_cofactors
            for p in reaction.left_participants + reaction.right_participants
        )

        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""

        # Exclude iron-sulfur cluster reactions (these are structural, not sulfur redox)
        iron_sulfur_patterns = ["2fe-2s", "4fe-4s", "[2fe-2s]", "[4fe-4s]",
                                "ferredoxin", "iron-sulfur"]
        if any(isp in label_lower for isp in iron_sulfur_patterns):
            return ClassificationResult(
                is_member=False,
                explanation="Iron-sulfur cluster reaction - not sulfur group redox"
            )

        # Key patterns for sulfur oxidoreductases (EC 1.8)
        # Must have explicit disulfide/dithiol interconversion or sulfite/sulfate redox

        # Check for explicit disulfide/dithiol patterns (both must be present)
        has_disulfide = "disulfide" in label_lower
        has_dithiol = "dithiol" in label_lower

        # True disulfide exchange: both disulfide and dithiol forms present
        is_disulfide_exchange = has_disulfide and has_dithiol

        # Glutathione reductase: requires BOTH glutathione disulfide AND reduced glutathione
        # Plus NAD(P)H cofactor - use ChEBI IDs as primary
        has_gssg = any(
            p.chebi_id == CHEBI_GLUTATHIONE_DISULFIDE
            for p in reaction.left_participants + reaction.right_participants
        ) or "glutathione disulfide" in label_lower

        has_gsh = any(
            p.chebi_id == CHEBI_GLUTATHIONE
            for p in reaction.left_participants + reaction.right_participants
        ) or ("glutathione" in label_lower and "disulfide" not in label_lower
              and "s-" not in label_lower)

        # Thioredoxin reductase: requires both disulfide and dithiol forms (check label)
        has_thioredoxin_ox = "thioredoxin" in label_lower and "disulfide" in label_lower
        has_thioredoxin_red = "thioredoxin" in label_lower and "dithiol" in label_lower

        # Glutathione reductase with NAD(P)H
        if has_gssg and has_gsh and has_redox:
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur oxidoreductase: glutathione reductase"
            )

        # Thioredoxin reductase
        if (has_thioredoxin_ox or has_thioredoxin_red) and has_redox:
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur oxidoreductase: thioredoxin reductase"
            )

        # Disulfide/dithiol exchange (even without NAD(P)H)
        if is_disulfide_exchange:
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur oxidoreductase: disulfide/dithiol exchange"
            )

        # Get label parts for substrate/product checking
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Sulfite/sulfate interconversion - but exclude monooxygenases that produce sulfite
        has_sulfite = any(
            p.chebi_id == CHEBI_SULFITE
            for p in reaction.left_participants + reaction.right_participants
        ) or "sulfite" in label_lower

        has_sulfate = any(
            p.chebi_id == CHEBI_SULFATE
            for p in reaction.left_participants + reaction.right_participants
        ) or ("sulfate" in label_lower and "thiosulfate" not in label_lower)

        # Exclude monooxygenases (O2 + sulfonate → aldehyde + sulfite)
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        ) or any(x in substrate_str for x in ["dioxygen", "oxygen", " o2"])
        has_aldehyde = "aldehyde" in product_str

        if has_o2 and has_aldehyde and has_sulfite:
            return ClassificationResult(
                is_member=False,
                explanation="Alkanesulfonate monooxygenase - not sulfur group redox"
            )

        # Sulfite/sulfate with redox cofactor indicates sulfur oxidoreductase
        if (has_sulfite or has_sulfate) and has_redox:
            return ClassificationResult(
                is_member=True,
                explanation="Sulfur oxidoreductase: sulfite/sulfate redox"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No sulfur oxidoreductase pattern"
        )
