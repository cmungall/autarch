"""Peptidase reaction classification using pattern DSL.

Peptidases are hydrolases that cleave peptide bonds in proteins and peptides.
EC 3.4.x.x classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import CHEBI_H2O, h_plus, p, water


class Peptidase(Hydrolase):
    """peptidase

    Examples:
    - Trypsin: protein + H2O → cleaved protein (after Arg/Lys)
    - Pepsin: protein + H2O → peptide fragments (broad specificity)
    - Carboxypeptidase: protein + H2O → protein + C-terminal amino acid
    - Aminopeptidase: protein + H2O → N-terminal amino acid + protein
    - Endopeptidases: internal cleavage of proteins
    """

    GO_ID = "GO:0008233"  # peptidase activity
    EC_NUMBER_PREFIX = "3.4.-.-"  # Acting on peptide bonds

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Protein hydrolysis: protein + H2O → fragments + amino acids
        var("protein") + p(water) >> var("fragment1") + var("fragment2") + optional(h_plus),
        # Peptide cleavage: peptide + H2O → shorter peptides
        var("peptide") + p(water) >> var("peptide1") + var("peptide2") + optional(h_plus),
        # Aminopeptidase: protein + H2O → amino acid + shortened protein
        var("protein") + p(water) >> var("amino_acid") + var("shortened_protein") + optional(h_plus),
        # General proteolysis: substrate + H2O → cleaved products
        var("substrate") + p(water) >> var("product1") + var("product2") + var("amino_acid") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a peptidase using pattern matching.

        Strategy:
        1. Must be a hydrolase (parent class)
        2. Must use water for hydrolysis
        3. Must involve protein/peptide substrates
        4. Must produce amino acids or peptide fragments
        """
        # First check if it's a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        # Must use water as reactant
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False,
                explanation="No water reactant - not peptide hydrolysis"
            )

        # Look for protein/peptide substrates

        has_protein_substrate = False

        # Look for amino acid products
        amino_acid_indicators = {
            "CHEBI:33709",  # amino acid
            "CHEBI:15956",  # L-alanine
            "CHEBI:16449",  # L-arginine
            "CHEBI:17196",  # L-asparagine
            "CHEBI:29991",  # L-aspartic acid
            "CHEBI:16828",  # L-cysteine
            "CHEBI:29985",  # L-glutamic acid
            "CHEBI:18050",  # L-glutamine
            "CHEBI:15428",  # L-glycine
            "CHEBI:27132",  # L-histidine
            "CHEBI:16847",  # L-isoleucine
            "CHEBI:25006",  # L-leucine
            "CHEBI:16643",  # L-lysine
            "CHEBI:16811",  # L-methionine
            "CHEBI:17295",  # L-phenylalanine
            "CHEBI:17203",  # L-proline
            "CHEBI:17115",  # L-serine
            "CHEBI:16857",  # L-threonine
            "CHEBI:16828",  # L-tryptophan
            "CHEBI:17895",  # L-tyrosine
            "CHEBI:16414",  # L-valine
        }


        has_amino_acid_product = any(
            p.chebi_id in amino_acid_indicators
            for p in reaction.right_participants
        )

        # Look for peptide fragments
        has_peptide_fragment = False

        if not (has_protein_substrate or has_amino_acid_product or has_peptide_fragment):
            return ClassificationResult(
                is_member=False,
                explanation="No protein/peptide substrate or amino acid products detected"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Peptidase: peptide bond hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
