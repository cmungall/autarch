"""Racemase/epimerase reaction classification using pattern DSL.

Racemases and epimerases catalyze stereoisomerization reactions.
EC 5.1.x.x classification.

This classifier inherits from Isomerase which handles:
- Water exclusion (hydrolases/lyases)
- Redox cofactor exclusion (NAD/NADP/FAD/O2)
- ATP/ADP balance checking
- CO2 exclusion (decarboxylases)
- Molecule count balance
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.stereochemistry import are_stereoisomers
from autarch.molecules import h_plus


class RacemaseAndEpimerase(Isomerase):
    """racemase and epimerase

    Examples:
    - Alanine racemase: L-alanine ⇌ D-alanine
    - Mandelate racemase: (S)-mandelate ⇌ (R)-mandelate
    - UDP-glucose 4-epimerase: UDP-glucose ⇌ UDP-galactose
    - Ribose-5-phosphate isomerase: ribose-5-P ⇌ ribulose-5-P
    """

    GO_ID = "GO:0016854"  # racemase and epimerase activity
    EC_NUMBER_PREFIX = "5.1.-.-"  # Racemases and epimerases

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple stereoisomerization: substrate ⇌ stereoisomer
        var("substrate") >> var("stereoisomer") + optional(h_plus),
        # Reversible reaction: L-form ⇌ D-form
        var("l_form") >> var("d_form"),
        # Sugar epimerase: UDP-sugar ⇌ UDP-epimer
        var("nucleotide_sugar") >> var("nucleotide_epimer") + optional(h_plus),
        # Amino acid racemase: L-amino acid ⇌ D-amino acid
        var("l_amino_acid") >> var("d_amino_acid") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a racemase/epimerase using pattern matching.

        Strategy:
        1. Must be an isomerase (parent class)
        2. Must be molecular formula conserving (no atoms added/removed)
        3. Must involve stereochemical change (L/D, R/S, alpha/beta)
        4. Look for characteristic substrate/product pairs
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # First check if it's an isomerase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an isomerase: {parent_result.explanation}",
            )

        diff = ReactionDiff(reaction)

        # Must be atom-conserving (stereoisomerization doesn't change formula)
        if not diff.is_balanced():
            return ClassificationResult(
                is_member=False,
                explanation="Not atom-balanced - not stereoisomerization"
            )

        # Must be 1:1 molecular conversion (substrate → product)
        if diff.n_reactant_molecules != 1 or diff.n_product_molecules != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Not 1:1 conversion - not simple stereoisomerization"
            )

        # Check for stereoisomerism using InChI (structural comparison)
        has_stereo_difference = False
        stereo_explanation = ""

        # Compare substrate and product using InChI for stereochemistry
        for reactant in reaction.left_participants:
            for product in reaction.right_participants:
                # Use InChI-based stereoisomer detection (preferred)
                if reactant.inchi and product.inchi:
                    is_stereo, explanation = are_stereoisomers(
                        reactant.inchi, product.inchi
                    )
                    if is_stereo:
                        has_stereo_difference = True
                        stereo_explanation = explanation
                        break

            if has_stereo_difference:
                break

        # Look for specific racemase/epimerase substrates by ChEBI ID
        racemase_chebi_pairs = {
            # L/D amino acid pairs
            ("CHEBI:16977", "CHEBI:15570"),  # L-alanine / D-alanine
            ("CHEBI:17822", "CHEBI:16344"),  # L-serine / D-serine
            ("CHEBI:17115", "CHEBI:16857"),  # L-glutamate / D-glutamate
            # UDP-sugar pairs
            ("CHEBI:18066", "CHEBI:18307"),  # UDP-glucose / UDP-galactose
        }

        # Check for known L/D or epimer pairs (one on each side)
        has_known_pair = False
        for l_id, d_id in racemase_chebi_pairs:
            left_ids = {p.chebi_id for p in reaction.left_participants}
            right_ids = {p.chebi_id for p in reaction.right_participants}
            if (l_id in left_ids and d_id in right_ids) or (d_id in left_ids and l_id in right_ids):
                has_known_pair = True
                break

        # Check for L/D or R/S naming patterns (stereoisomer indicators in names)
        stereo_name_patterns = [
            ("l-", "d-"),  # L-alanine / D-alanine
            ("-l-", "-d-"),  # alpha-L-arabinose / alpha-D-xylose
            ("(s)-", "(r)-"),  # (S)-mandelate / (R)-mandelate
            ("(s)", "(r)"),  # alternate parentheses format
            ("alpha-", "beta-"),  # alpha-D-glucose / beta-D-glucose
        ]

        # Sugar epimer pairs (glucose/galactose, glucosamine/galactosamine, etc.)
        # These are epimerizations at specific carbons, not simple L/D changes
        sugar_epimer_pairs = [
            ("glucose", "galactose"),  # C4 epimer
            ("glucuronate", "galacturonate"),  # C4 epimer
            ("glucosamine", "galactosamine"),  # C4 epimer
            ("glucosamine", "mannosamine"),  # C2 epimer
            ("glucose", "mannose"),  # C2 epimer
            ("xylose", "arabinose"),  # C4 epimer
        ]

        # Use label strings for pattern matching
        has_stereo_naming = False
        for left_prefix, right_prefix in stereo_name_patterns:
            # Check if one side has L-/S-/alpha- and other has D-/R-/beta-
            if (left_prefix in substrate_str and right_prefix in product_str) or \
               (right_prefix in substrate_str and left_prefix in product_str):
                has_stereo_naming = True
                break

        # Check for sugar epimer pairs
        has_epimer_pair = False
        for sugar1, sugar2 in sugar_epimer_pairs:
            if (sugar1 in substrate_str and sugar2 in product_str) or \
               (sugar2 in substrate_str and sugar1 in product_str):
                has_epimer_pair = True
                break

        # Require actual evidence of stereoisomerization
        if not (has_stereo_difference or has_known_pair or has_stereo_naming or has_epimer_pair):
            return ClassificationResult(
                is_member=False,
                explanation="No stereochemical conversion detected (no InChI stereo diff, no L/D pair, no stereo naming, no sugar epimer pair)"
            )

        # Apply racemase-specific exclusions
        # (Parent Isomerase already excludes: water, NAD/NADP/FAD, ATP imbalance, CO2)

        # Exclude aldose-ketose isomerases (constitutional, not stereoisomers)
        # E.g., glucose 6-P → fructose 6-P (ring-open aldose → ketose)
        aldose_ketose_pairs = [
            ("glucose", "fructose"),
            ("ribose", "ribulose"),
            ("xylose", "xylulose"),
            ("mannose", "fructose"),
        ]
        for aldose, ketose in aldose_ketose_pairs:
            if (aldose in substrate_str and ketose in product_str) or \
               (ketose in substrate_str and aldose in product_str):
                return ClassificationResult(
                    is_member=False,
                    explanation="Aldose-ketose isomerase - constitutional change, not stereoisomerization"
                )

        # Exclude phosphomutases (1-phosphate → 6-phosphate positional change)
        # E.g., GlcNAc 1-phosphate → GlcNAc 6-phosphate
        has_1_phosphate_left = "1-phosphate" in substrate_str
        has_6_phosphate_right = "6-phosphate" in product_str
        has_6_phosphate_left = "6-phosphate" in substrate_str
        has_1_phosphate_right = "1-phosphate" in product_str
        if (has_1_phosphate_left and has_6_phosphate_right) or (has_6_phosphate_left and has_1_phosphate_right):
            return ClassificationResult(
                is_member=False,
                explanation="Phosphomutase - positional isomer, not stereoisomer"
            )

        # Exclude double bond isomerases (different double bond positions)
        # E.g., lathosterol → cholest-8-en-ol (different double bond location in steroid)
        # Check if label has different -en- positions (double bond isomers)
        left_en_pos = re.findall(r'-(\d+)-en-', substrate_str)
        right_en_pos = re.findall(r'-(\d+)-en-', product_str)
        if left_en_pos and right_en_pos and left_en_pos != right_en_pos:
            return ClassificationResult(
                is_member=False,
                explanation="Double bond position isomerase - constitutional change, not stereoisomer"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Determine specific type based on substrates
        explanation = "Racemase/epimerase: stereoisomerization"

        # Check for amino acid racemase
        amino_acid_patterns = ["alanine", "serine", "glutamate", "aspartate", "proline", "methionine"]
        if any(aa in label_lower for aa in amino_acid_patterns):
            explanation = "Racemase/epimerase: amino acid racemization"

        # Check for UDP-sugar epimerase or sugar epimer pairs
        if any(prefix in label_lower for prefix in ["udp-", "gdp-", "dtdp-"]):
            explanation = "Racemase/epimerase: nucleotide-sugar epimerization"
        elif has_epimer_pair:
            explanation = "Racemase/epimerase: sugar epimerization"

        # Check for mandelate/lactate racemase
        if "mandelate" in label_lower:
            explanation = "Racemase/epimerase: mandelate racemization"
        elif "lactate" in label_lower:
            explanation = "Racemase/epimerase: lactate racemization"

        # Include stereochemistry details from InChI analysis
        if has_stereo_difference and stereo_explanation:
            explanation += f" [{stereo_explanation}]"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
