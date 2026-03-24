"""CH-NH Oxidoreductase reaction classification.

EC 1.5: Oxidoreductases acting on the CH-NH group of donors.
These enzymes catalyze oxidation of secondary amines (C-NH-C) or primary amines
attached to secondary carbons.

Pattern: R-NH-R' + acceptor → R=N-R' + reduced acceptor
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_H2O2,
    CHEBI_H2O,
)

# Common CH-NH oxidoreductase substrates
CHEBI_PROLINE = "CHEBI:26271"
CHEBI_SARCOSINE = "CHEBI:15611"
CHEBI_GLYCINE_BETAINE = "CHEBI:17750"
CHEBI_DIMETHYLGLYCINE = "CHEBI:17724"
CHEBI_TETRAHYDROFOLATE = "CHEBI:57453"

# Products characteristic of CH-NH oxidation
CHEBI_PYRROLINE_CARBOXYLATE = "CHEBI:17388"  # From proline
CHEBI_FORMALDEHYDE = "CHEBI:16842"


class OxidoreductaseActingOnTheCHNHGroupOfDonors(Oxidoreductase):
    """oxidoreductase acting on the CH-NH group of donors

    EC 1.5.x.x includes:
    - Proline dehydrogenase (proline → pyrroline-5-carboxylate)
    - Sarcosine dehydrogenase (sarcosine → glycine + formaldehyde)
    - Dimethylglycine dehydrogenase
    - Methylenetetrahydrofolate reductase
    - Dihydrofolate reductase (dihydrofolate → tetrahydrofolate)

    Examples:
    - L-proline + acceptor → 1-pyrroline-5-carboxylate + reduced acceptor
    - Sarcosine + FAD → glycine + formaldehyde + FADH2
    """

    GO_ID = "GO:0016645"  # oxidoreductase activity, acting on the CH-NH group
    EC_NUMBER_PREFIX = "1.5.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a CH-NH oxidoreductase.

        Strategy:
        1. Must be an oxidoreductase (have redox cofactor)
        2. Look for characteristic CH-NH substrates/products
        3. Exclude CH-NH2 oxidoreductases (EC 1.4) by product patterns
        """
        # Must have redox cofactor
        redox_cofactors = {
            CHEBI_NAD_PLUS, CHEBI_NADH, CHEBI_NADP_PLUS, CHEBI_NADPH,
            CHEBI_FAD, CHEBI_FADH2, CHEBI_O2, CHEBI_H2O2,
            "CHEBI:58210", "CHEBI:58307",  # FMN/FMNH2
        }

        has_redox = any(
            p.chebi_id in redox_cofactors
            for p in reaction.left_participants + reaction.right_participants
        )

        if not has_redox:
            return ClassificationResult(
                is_member=False,
                explanation="No redox cofactor - not oxidoreductase"
            )

        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""

        # Get substrates (left side of equation) from label
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Look for characteristic CH-NH substrates
        # These are secondary amines or compounds with C-NH-C bonds
        chnh_substrate_patterns = [
            "proline", "pyrroline", "pyrrolidine", "piperidine",
            "sarcosine", "dimethylglycine", "betaine",
            "dihydrofolate", "tetrahydrofolate", "methylenetetrahydro",
            "formyltetrahydro", "methenyltetrahydro",
            "trimethylamine", "dimethylamine", "methylamine",
            "saccharopin", "pipecol",
            "pyrimidine", "imidazole", "histidine",
            "spermidine", "spermine", "putrescine",
            "nicotinamide",
            "nicotine", "nicotinium", "nornicotine", "anatabine",
            "-hydroxynicotine", "-hydroxynicotinium",
            "berberine", "reticuline",
            "opine", "octopine", "tauropine", "strombine",
            "flavin",
        ]

        has_chnh_substrate = any(pat in substrate_str for pat in chnh_substrate_patterns)

        # Check for characteristic products (use label, not p.name)
        chnh_product_patterns = [
            "pyrroline", "imine",
            "formaldehyde", "formate",
            "tetrahydrofolate",
            "ammonia", "nh4", "nh3",
        ]

        any(pat in product_str for pat in chnh_product_patterns)

        # Check for formaldehyde production (common in CH-NH oxidation of N-methyl compounds)
        # Primary: ChEBI ID; Fallback: label pattern
        has_formaldehyde = any(
            p.chebi_id == CHEBI_FORMALDEHYDE
            for p in reaction.right_participants
        ) or "formaldehyde" in product_str

        # Exclude simple amino acid dehydrogenases (EC 1.4 - act on CH-NH2)
        # These typically have alpha-keto acid products from amino acid substrates
        keto_acid_patterns = ["2-oxo", "alpha-keto", "α-keto", "oxaloacet"]
        any(pat in product_str for pat in keto_acid_patterns)

        # EC 1.4 amino acid oxidoreductases typically convert R-CH(NH2)-COOH to R-C(=O)-COOH
        # We want to exclude these (check in full label)
        amino_acid_dehydrogenase_patterns = [
            "glutamate dehydrogen", "alanine dehydrogen", "aspartate dehydrogen",
            "leucine dehydrogen", "valine dehydrogen", "amino acid dehydrogen"
        ]

        if any(pat in label_lower for pat in amino_acid_dehydrogenase_patterns):
            return ClassificationResult(
                is_member=False,
                explanation="Amino acid dehydrogenase (EC 1.4) - not CH-NH group"
            )

        # Exclude transaminases/aminotransferases (EC 2.6) - they transfer NH2, not oxidize CH-NH
        if "transaminase" in label_lower:
            return ClassificationResult(
                is_member=False,
                explanation="Transaminase - not CH-NH oxidoreductase"
            )

        # Exclude monooxygenases (EC 1.14) - they hydroxylate using O2
        # Pattern: substrate + O2 + reductant → hydroxylated-substrate + H2O
        # Primary: ChEBI ID; Fallback: label pattern
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        ) or any(x in substrate_str for x in ["dioxygen", "oxygen", " o2"])

        # Monooxygenases produce H2O alongside hydroxylated product
        produces_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.right_participants
        )

        # Check for hydroxylated product pattern (in label)
        has_hydroxylated_product = "hydroxy" in product_str or "-oh" in product_str

        # Monooxygenase pattern: O2 consumed, water produced, product hydroxylated
        if has_o2 and produces_water and has_hydroxylated_product:
            return ClassificationResult(
                is_member=False,
                explanation="Monooxygenase (EC 1.14) - not CH-NH oxidoreductase"
            )

        # Match if we have characteristic CH-NH patterns
        if has_chnh_substrate:
            # Find which pattern matched for explanation
            matched_pattern = next(
                (pat for pat in chnh_substrate_patterns if pat in substrate_str),
                "secondary amine"
            )
            return ClassificationResult(
                is_member=True,
                explanation=f"CH-NH oxidoreductase: {matched_pattern} oxidation"
            )

        # Formaldehyde + specific N-methyl substrate patterns
        # (not formaldehyde alone, as that can come from S-chemistry)
        if has_formaldehyde and has_redox:
            n_methyl_patterns = ["methylamine", "dimethyl", "trimethyl", "n-methyl"]
            if any(pat in substrate_str for pat in n_methyl_patterns):
                return ClassificationResult(
                    is_member=True,
                    explanation="CH-NH oxidoreductase: N-demethylation"
                )

        return ClassificationResult(
            is_member=False,
            explanation="No CH-NH oxidoreductase pattern"
        )
