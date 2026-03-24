"""Phosphoric diester hydrolase reaction classification.

Phosphoric diester hydrolases (EC 3.1.4) are hydrolases acting on phosphoric diester bonds.
They cleave P-O-C bonds in phosphodiesters.

Pattern: phosphodiester + H2O → alcohol + phosphomonoester
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.molecules import (
    CHEBI_H2O,
)


class PhosphoricDiesterHydrolase(Hydrolase):
    """Classifier for phosphoric diester hydrolase reactions.

    EC 3.1.4.x includes:
    - Phospholipase C: phosphatidylcholine → diacylglycerol + phosphocholine
    - Phospholipase D: phosphatidylcholine → phosphatidate + choline
    - Cyclic nucleotide phosphodiesterase: cAMP → AMP
    - Sphingomyelin phosphodiesterase

    Examples:
    - Phosphatidylcholine + H2O → DAG + phosphocholine
    - cAMP + H2O → AMP
    """

    GO_ID = "GO:0008081"  # phosphoric diester hydrolase activity
    EC_NUMBER_PREFIX = "3.1.4.-"  # Phosphoric diester hydrolases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a phosphoric diester hydrolase.

        Strategy:
        1. Must use water as substrate
        2. Look for phosphodiester substrates (phospholipids, cyclic nucleotides)
        3. Products include phosphomonoester + alcohol OR cyclic → linear
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must use water
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water substrate - phosphoric diester hydrolases are hydrolases"
            )

        # Check for phosphodiester substrates by label patterns
        phosphodiester_patterns = [
            "phosphatidyl", "phospholipid", "sphingomyelin",
            "cyclic amp", "cyclic gmp", "camp", "cgmp",
            "glycero-3-phospho", "phosphocholine",
            "diacyl-sn-glycero-3-phospho"
        ]

        has_phosphodiester = any(pat in substrate_str for pat in phosphodiester_patterns)

        # Check for phospholipase products
        phospholipase_products = [
            "diacylglycerol", "diglyceride", "dag",
            "phosphatidate", "phosphatidic acid",
            "phosphocholine", "phosphoethanolamine",
            "choline", "ethanolamine",
            "ceramide", "sphingosine"
        ]

        has_phospholipase_product = any(pat in product_str for pat in phospholipase_products)

        # Check for cyclic nucleotide phosphodiesterase pattern
        # cAMP → AMP, cGMP → GMP
        cyclic_nucleotides = ["cyclic amp", "cyclic gmp", "camp", "cgmp",
                             "3',5'-cyclic", "adenosine 3',5'-cyclic"]
        linear_nucleotides = ["amp", "gmp", "adenosine 5'-monophosphate",
                             "guanosine 5'-monophosphate"]

        has_cyclic_substrate = any(cn in substrate_str for cn in cyclic_nucleotides)

        has_linear_product = any(ln in product_str for ln in linear_nucleotides)

        if has_phosphodiester and has_phospholipase_product:
            # Determine specific phospholipase type
            if "diacylglycerol" in product_str:
                return ClassificationResult(
                    is_member=True,
                    explanation="Phosphoric diester hydrolase: phospholipase C-type (-> DAG + phospho-alcohol)"
                )
            elif "phosphatidat" in product_str:
                return ClassificationResult(
                    is_member=True,
                    explanation="Phosphoric diester hydrolase: phospholipase D-type (-> phosphatidate + alcohol)"
                )
            else:
                return ClassificationResult(
                    is_member=True,
                    explanation="Phosphoric diester hydrolase: phospholipid hydrolysis"
                )

        if has_cyclic_substrate and has_linear_product:
            return ClassificationResult(
                is_member=True,
                explanation="Phosphoric diester hydrolase: cyclic nucleotide -> linear nucleotide"
            )

        # Generic phosphodiester pattern
        if has_phosphodiester:
            return ClassificationResult(
                is_member=True,
                explanation="Phosphoric diester hydrolase: phosphoric diester hydrolysis"
            )

        # Check for glycerophosphodiester substrates (NOT monoesters like "glycerol 2-phosphate")
        # Phosphodiesters have patterns like "glycero-3-phospho-X" or "phosphodiester"
        glycerophosphodiester_patterns = [
            "glycero-3-phospho",  # sn-glycero-3-phosphocholine, etc.
            "phosphodiester",
            "glycerophospho",  # glycerophosphocholine, etc.
        ]

        has_glycerophosphodiester = any(pat in substrate_str for pat in glycerophosphodiester_patterns)

        if has_glycerophosphodiester:
            return ClassificationResult(
                is_member=True,
                explanation="Phosphoric diester hydrolase: glycerophosphodiester hydrolysis"
            )

        # Note: "glycerol phosphate" alone is a MONOester (phosphatase substrate),
        # not a DIester. Don't match simple glycerol + phosphate patterns.

        return ClassificationResult(
            is_member=False,
            explanation="No phosphodiester hydrolysis pattern"
        )
