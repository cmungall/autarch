"""Terpene synthase reaction classification.

Terpene synthases (EC 4.2.3) are carbon-oxygen lyases acting on phosphates.
They catalyze the formation of terpenes/terpenoids from prenyl diphosphate precursors
by releasing diphosphate (pyrophosphate).

Pattern: prenyl-diphosphate → terpene/terpenoid + diphosphate
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.lyase import Lyase
from autarch.molecules import (
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
)

# Prenyl diphosphate precursors
CHEBI_DMAPP = "CHEBI:57623"  # dimethylallyl diphosphate
CHEBI_IPP = "CHEBI:58536"  # isopentenyl diphosphate
CHEBI_GPP = "CHEBI:58057"  # geranyl diphosphate
CHEBI_FPP = "CHEBI:175763"  # farnesyl diphosphate (2-trans,6-trans)
CHEBI_GGPP = "CHEBI:57533"  # geranylgeranyl diphosphate

# Set of all prenyl diphosphate precursors
PRENYL_DIPHOSPHATES = {
    CHEBI_DMAPP,
    CHEBI_IPP,
    CHEBI_GPP,
    CHEBI_FPP,
    CHEBI_GGPP,
    "CHEBI:58756",  # (2E,6E,10E)-geranylgeranyl diphosphate
    "CHEBI:57533",  # GGPP
}


class TerpeneSynthase(Lyase):
    """terpene synthase

    These enzymes (EC 4.2.3.x) catalyze C-O bond cleavage to release diphosphate
    while forming carbocyclic or linear terpene products.

    Examples:
    - Limonene synthase: geranyl-PP → limonene + PPi
    - Vetispiradiene synthase: farnesyl-PP → vetispiradiene + PPi
    - Copalyl diphosphate synthase: GGPP → copalyl-PP
    """

    GO_ID = "GO:0010333"  # terpene synthase activity
    EC_NUMBER_PREFIX = "4.2.3.-"  # Carbon-oxygen lyases acting on phosphates
    EC_BROAD_XREFS = ["4.2.3.-"]  # GO xref is broadMatch, not exact EC equivalence

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a terpene synthase.

        Strategy:
        1. Must produce diphosphate as a product
        2. Must have prenyl diphosphate precursor as substrate
        3. Must be single substrate (cyclization) - not prenyltransferase
        4. Product must NOT be a prenyl-PP (not chain elongation)
        5. Exclude hydrolases (water substrates)
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Check for diphosphate product (hallmark of terpene synthases)
        has_ppi_product = any(
            p.chebi_id == CHEBI_DIPHOSPHATE
            for p in reaction.right_participants
        )

        if not has_ppi_product:
            return ClassificationResult(
                is_member=False,
                explanation="No diphosphate product - terpene synthases release PPi"
            )

        # Check for prenyl diphosphate substrate - ChEBI ID primary
        prenyl_substrates = [
            p for p in reaction.left_participants
            if p.chebi_id in PRENYL_DIPHOSPHATES
        ]

        # Fallback: check label for prenyl patterns
        prenyl_patterns = ["geranyl", "farnesyl", "isopentenyl", "dimethylallyl"]
        has_prenyl_in_label = (
            "diphosphate" in substrate_str and
            any(pn in substrate_str for pn in prenyl_patterns)
        )

        if not prenyl_substrates and not has_prenyl_in_label:
            return ClassificationResult(
                is_member=False,
                explanation="No prenyl diphosphate substrate"
            )

        # Exclude hydrolases (use water as substrate)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Uses water as substrate - hydrolase/phosphatase"
            )

        # Exclude chain elongation (IPP + DMAPP → GPP + PPi)
        # These produce prenyl-diphosphate products - check ChEBI ID primary
        has_prenyl_product = any(
            p.chebi_id in PRENYL_DIPHOSPHATES
            for p in reaction.right_participants
        )

        # Fallback: check label for prenyl product patterns
        if not has_prenyl_product:
            prenyl_product_patterns = ["geranyl", "farnesyl", "neryl"]
            has_prenyl_product = (
                "diphosphate" in product_str and
                any(pn in product_str for pn in prenyl_product_patterns)
            )

        if has_prenyl_product:
            return ClassificationResult(
                is_member=False,
                explanation="Produces prenyl-diphosphate - chain elongation prenyltransferase"
            )

        # Count non-prenyl substrates (acceptors for prenyltransferases)
        non_prenyl_substrates = [
            p for p in reaction.left_participants
            if p not in prenyl_substrates and
            p.chebi_id != CHEBI_H2O and
            p.chebi_id != CHEBI_DIPHOSPHATE
        ]

        # Terpene synthases should have mainly prenyl-PP substrates
        # Prenyltransferases have an acceptor molecule (protein, amino acid, nucleotide, etc.)
        if len(non_prenyl_substrates) >= 1:
            # Known acceptor ChEBI IDs for prenyltransferases
            # (Many compound names are truncated SMILES, so use ChEBI IDs as primary detection)
            ACCEPTOR_CHEBI = {
                # Amino acids (common prenylation targets)
                "CHEBI:57912",  # L-tryptophan
                "CHEBI:58315",  # L-tyrosine
                "CHEBI:57972",  # L-phenylalanine
                "CHEBI:57595",  # L-histidine
                "CHEBI:58359",  # L-cysteine
                # Aromatic compounds
                "CHEBI:17879",  # 4-hydroxybenzoate
                "CHEBI:36059",  # benzoate
                "CHEBI:28966",  # 1,4-dihydroxynaphthalene
                # Nucleosides/nucleotides
                "CHEBI:16335",  # adenosine
                "CHEBI:58245",  # AMP
                "CHEBI:456215",  # ADP
                # tRNA/polymer indicators
                "CHEBI:74411",  # generic polymer (often tRNA)
            }

            has_acceptor_chebi = any(
                p.chebi_id in ACCEPTOR_CHEBI
                for p in non_prenyl_substrates
            )

            # Acceptor patterns for label-based detection (fallback when ChEBI not in set)
            acceptor_patterns = [
                # Proteins and peptides
                "protein", "peptide", "cysteinyl", "cysteine", "trna", "rrna",
                "polymer",
                # Amino acids
                "tryptophan", "tyrosine", "phenylalanine", "histidine",
                "alanine", "glycine", "serine", "threonine",
                "leucine", "isoleucine", "valine", "proline",
                "methionine", "asparagine", "glutamine",
                "aspartate", "glutamate", "lysine", "arginine",
                # Nucleotides and nucleosides
                "adenosine", "guanosine", "cytidine", "uridine", "thymidine",
                "adenine", "guanine", "cytosine", "uracil", "thymine",
                # Lipids and glycerol
                "glycerol", "glyceride", "lipid", "phospholipid",
                # Aromatic compounds
                "pterocarp", "flavon", "coumarin", "phenyl", "benzyl",
                "benzoate", "hydroxybenzoate", "hydroxybenzoic",
                "indole", "quinone", "naphtho", "naphthol",
                # Sugars
                "glucose", "galactose", "mannose", "ribose",
            ]

            # Check label for acceptor patterns (not p.name)
            has_acceptor_label = any(pat in substrate_str for pat in acceptor_patterns)

            if has_acceptor_chebi or has_acceptor_label:
                # Find which acceptor pattern matched
                matched_acceptor = next(
                    (pat for pat in acceptor_patterns if pat in substrate_str),
                    "unknown"
                )
                return ClassificationResult(
                    is_member=False,
                    explanation=f"Has acceptor substrate ({matched_acceptor}) - prenyltransferase not terpene synthase"
                )

        # Multiple prenyl substrates suggests chain elongation (IPP + GPP → FPP)
        if len(prenyl_substrates) >= 2:
            return ClassificationResult(
                is_member=False,
                explanation="Multiple prenyl substrates - chain elongation prenyltransferase"
            )

        # Determine terpene type based on substrate - ChEBI ID primary, label fallback
        terpene_type = "unknown"
        for participant in reaction.left_participants:
            if participant.chebi_id == CHEBI_GPP:
                terpene_type = "monoterpene (C10)"
                break
            elif participant.chebi_id == CHEBI_FPP:
                terpene_type = "sesquiterpene (C15)"
                break
            elif participant.chebi_id == CHEBI_GGPP:
                terpene_type = "diterpene (C20)"
                break
            elif participant.chebi_id in {CHEBI_DMAPP, CHEBI_IPP}:
                terpene_type = "hemiterpene (C5)"
                break

        # Fallback to label patterns
        if terpene_type == "unknown":
            if "geranyl" in substrate_str and "geranylgeranyl" not in substrate_str:
                terpene_type = "monoterpene (C10)"
            elif "farnesyl" in substrate_str:
                terpene_type = "sesquiterpene (C15)"
            elif "geranylgeranyl" in substrate_str:
                terpene_type = "diterpene (C20)"

        explanation = f"Terpene synthase: prenyl-PP → {terpene_type} + PPi"

        return ClassificationResult(
            is_member=True,
            explanation=explanation
        )
