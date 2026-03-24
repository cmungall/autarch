"""Prenyltransferase reaction classification.

Prenyltransferases (EC 2.5.1) transfer alkyl or aryl groups other than methyl.
The most common are prenyltransferases that transfer prenyl groups (C5 units)
from prenyl diphosphate donors.

Pattern: prenyl-PP + acceptor → prenylated-product + PPi
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.transferase import Transferase
from autarch.molecules import (
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_ATP,
    CHEBI_ADP,
    CHEBI_AMP,
)

# Prenyl diphosphate donors
CHEBI_DMAPP = "CHEBI:57623"  # dimethylallyl diphosphate
CHEBI_IPP = "CHEBI:58536"  # isopentenyl diphosphate
CHEBI_GPP = "CHEBI:58057"  # geranyl diphosphate
CHEBI_FPP = "CHEBI:175763"  # farnesyl diphosphate (2-trans,6-trans)
CHEBI_GGPP = "CHEBI:57533"  # geranylgeranyl diphosphate

PRENYL_DIPHOSPHATES = {
    CHEBI_DMAPP,
    CHEBI_IPP,
    CHEBI_GPP,
    CHEBI_FPP,
    CHEBI_GGPP,
    "CHEBI:58756",  # (2E,6E,10E)-geranylgeranyl diphosphate
}


class Prenyltransferase(Transferase):
    """prenyltransferase

    Most commonly, these are prenyl transfers from isoprenoid diphosphates.
    EC 2.5.1.x includes:
    - Farnesyltransferase
    - Geranylgeranyltransferase
    - Dimethylallyltransferase

    Examples:
    - DMAPP + acceptor → prenylated-acceptor + PPi
    - FPP + protein → farnesylated-protein + PPi
    """

    GO_ID = "GO:0004659"  # prenyltransferase activity
    EC_NUMBER_PREFIX = "2.5.1.-"  # Transferring alkyl or aryl groups
    EC_BROAD_XREFS = ["2.5.1.-"]  # GO xref is broadMatch, not exact EC equivalence

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a prenyltransferase.

        Strategy:
        1. Look for prenyl diphosphate donor on left
        2. Look for diphosphate product on right
        3. Exclude hydrolases (water as substrate)
        4. Exclude terpene synthases (single prenyl substrate → cyclic terpene)
        5. Allow chain elongation (IPP + prenyl-PP → longer prenyl-PP + PPi)
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        prenyl_names = ["dimethylallyl", "isopentenyl", "geranyl", "farnesyl", "geranylgeranyl"]

        # Identify prenyl diphosphate substrates - ChEBI ID primary
        prenyl_substrates = [
            p for p in reaction.left_participants
            if p.chebi_id in PRENYL_DIPHOSPHATES
        ]

        # Fallback: check label for prenyl patterns
        has_prenyl_in_label = (
            "diphosphate" in substrate_str and
            any(pn in substrate_str for pn in prenyl_names)
        )

        if not prenyl_substrates and not has_prenyl_in_label:
            return ClassificationResult(
                is_member=False,
                explanation="No prenyl diphosphate donor substrate"
            )

        # Check for diphosphate product
        has_ppi_product = any(
            p.chebi_id == CHEBI_DIPHOSPHATE
            for p in reaction.right_participants
        )

        if not has_ppi_product:
            return ClassificationResult(
                is_member=False,
                explanation="No diphosphate product - prenyltransferases release PPi"
            )

        # Exclude hydrolases (water as substrate converts prenyl-PP to prenyl alcohol)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Uses water - phosphatase/hydrolase not prenyltransferase"
            )

        # Exclude ATP-dependent phosphotransferases
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        has_adp = any(p.chebi_id in {CHEBI_ADP, CHEBI_AMP} for p in reaction.right_participants)

        if has_atp and has_adp:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphotransferase pattern - ATP dependent"
            )

        # Identify non-prenyl substrates (potential acceptors)
        non_prenyl_substrates = [
            p for p in reaction.left_participants
            if p not in prenyl_substrates and
            p.chebi_id != CHEBI_H2O
        ]

        # Check if this is chain elongation (two prenyl substrates → longer prenyl-PP + PPi)
        # IPP + DMAPP → GPP + PPi, IPP + GPP → FPP + PPi, etc.
        if len(prenyl_substrates) >= 2 or (prenyl_substrates and has_prenyl_in_label):
            # Check for prenyl-PP product (chain elongation) - ChEBI ID primary
            has_prenyl_product = any(
                p.chebi_id in PRENYL_DIPHOSPHATES
                for p in reaction.right_participants
            )
            # Fallback to label
            if not has_prenyl_product:
                has_prenyl_product = (
                    "diphosphate" in product_str and
                    any(pn in product_str for pn in prenyl_names)
                )
            if has_prenyl_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="Prenyltransferase: prenyl chain elongation"
                )

        # Exclude terpene synthases (only prenyl substrates → cyclic product)
        # Terpene synthases have only 1 prenyl substrate and no acceptor
        if len(prenyl_substrates) == 1 and len(non_prenyl_substrates) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Terpene synthase pattern - single prenyl substrate, no acceptor"
            )

        # Must have an acceptor to be a prenyltransferase
        if len(non_prenyl_substrates) == 0 and len(prenyl_substrates) == 1:
            return ClassificationResult(
                is_member=False,
                explanation="No acceptor substrate detected"
            )

        # Determine prenyl type - ChEBI ID primary, label fallback
        prenyl_type = "alkyl/aryl"
        for participant in reaction.left_participants:
            if participant.chebi_id == CHEBI_DMAPP:
                prenyl_type = "dimethylallyl (C5)"
                break
            elif participant.chebi_id == CHEBI_GPP:
                prenyl_type = "geranyl (C10)"
                break
            elif participant.chebi_id == CHEBI_FPP:
                prenyl_type = "farnesyl (C15)"
                break
            elif participant.chebi_id == CHEBI_GGPP:
                prenyl_type = "geranylgeranyl (C20)"
                break

        # Fallback to label patterns
        if prenyl_type == "alkyl/aryl":
            if "dimethylallyl" in substrate_str:
                prenyl_type = "dimethylallyl (C5)"
            elif "geranyl" in substrate_str and "geranylgeranyl" not in substrate_str:
                prenyl_type = "geranyl (C10)"
            elif "farnesyl" in substrate_str:
                prenyl_type = "farnesyl (C15)"
            elif "geranylgeranyl" in substrate_str:
                prenyl_type = "geranylgeranyl (C20)"

        return ClassificationResult(
            is_member=True,
            explanation=f"Prenyltransferase: {prenyl_type} group transfer"
        )
