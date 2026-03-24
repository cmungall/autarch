"""Amidase reaction classification using pattern DSL.

Amidases are hydrolases that cleave amide bonds.
EC 3.5.x.x classification.

History
-------

## 2025-12-22

Expanded nitrogen_products set with additional ChEBI IDs from RHEA:
- CHEBI:28938 (ammonium) - 132 occurrences in RHEA
- More amino acid product ChEBI IDs
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_H2O,
    CHEBI_NH3,
    CHEBI_NH4,
    h_plus,
    p,
    water,
)


class HydrolaseActingOnCarbonNitrogenButNotPeptideBonds(Hydrolase):
    """hydrolase acting on carbon-nitrogen but not peptide bonds

    Examples:
    - Acetamidase: acetamide + H2O → acetate + NH3
    - Formamidase: formamide + H2O → formate + NH3
    - Fatty acid amide hydrolase: fatty acid amide + H2O → fatty acid + amine
    - Glutaminase: glutamine + H2O → glutamate + NH3
    - Asparagine: asparagine + H2O → aspartate + NH3
    """

    GO_ID = "GO:0016810"  # hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds
    EC_NUMBER_PREFIX = "3.5.-.-"  # Acting on C-N bonds, other than peptide bonds
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")  # generic amide C(=O)-N
    CHEBI_SULFATE = "CHEBI:16189"

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple amide hydrolysis: amide + H2O → carboxylate + amine/NH3
        var("amide") + p(water) >> var("carboxylate") + var("amine") + optional(h_plus),
        # Amino acid deamidation: amino acid amide + H2O → amino acid + NH3
        var("amino_acid_amide") + p(water) >> var("amino_acid") + var("ammonia") + optional(h_plus),
        # Fatty acid amide: fatty acid amide + H2O → fatty acid + amine
        var("fatty_amide") + p(water) >> var("fatty_acid") + var("amine") + optional(h_plus),
        # Simple amide: acetamide + H2O → acetate + NH3 + H+
        var("simple_amide") + p(water) >> var("acid") + var("ammonia") + p(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an amidase using pattern matching.

        Strategy:
        1. Must be a hydrolase (parent class)
        2. Must use water for hydrolysis
        3. Must involve amide bond cleavage (C-N bond)
        4. Must produce amine/ammonia as product
        5. Exclude peptidases (which are specialized amidases)
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
                explanation="No water reactant - not amide hydrolysis"
            )

        # Exclude thioester hydrolysis (CoA/thiol esters): these are C-S hydrolases.
        has_coa = any(
            p.chebi_id == CHEBI_COA
            for p in reaction.left_participants + reaction.right_participants
        )
        has_thioester_substrate = any(p.is_thioester() for p in reaction.left_participants)
        if has_coa or has_thioester_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="CoA/thioester hydrolysis - C-S bond chemistry, not amidase"
            )

        # Exclude sulfatase chemistry (sulfate release from sulfate esters).
        has_sulfate_product = any(
            p.chebi_id == self.CHEBI_SULFATE
            for p in reaction.right_participants
        )
        if has_sulfate_product:
            return ClassificationResult(
                is_member=False,
                explanation="Sulfate release indicates sulfatase chemistry, not amidase"
            )

        # Look for amide substrates or amide-derived products
        amide_indicators = {
            # Common amides
            "CHEBI:27856",  # acetamide
            "CHEBI:16327",  # formamide
            "CHEBI:18050",  # glutamine
            "CHEBI:17196",  # asparagine
            "CHEBI:18307",  # carbamoyl phosphate (contains amide)
        }

        # Check for amide substrates by ChEBI ID
        has_amide_substrate = any(
            p.chebi_id in amide_indicators
            for p in reaction.left_participants
        )
        has_amide_substrate = has_amide_substrate or any(
            self._has_amide_moiety(p) for p in reaction.left_participants
        )

        # Look for products indicating amide bond cleavage
        nitrogen_products = {
            CHEBI_NH3,      # ammonia (CHEBI:16134)
            CHEBI_NH4,      # ammonium (CHEBI:28938)
            # Amino acid products from deamidation
            "CHEBI:29985",  # L-glutamate (from L-glutamine)
            "CHEBI:29991",  # L-aspartate (from L-asparagine)
            "CHEBI:16015",  # L-glutamate(2-)
            "CHEBI:29958",  # L-aspartate(1-)
            # Common amine products
            "CHEBI:15571",  # L-cysteine (from asparagine-linked)
            "CHEBI:16467",  # L-arginine
            "CHEBI:58048",  # L-ornithine
        }

        # Products indicating N-deacetylation (acetate released)
        # CRITICAL: Acetate alone is NOT sufficient - it could be ESTER hydrolysis!
        # Amidase: R-NH-CO-CH3 + H2O → R-NH2 + acetate (produces amine + acetate)
        # Esterase: R-O-CO-CH3 + H2O → R-OH + acetate (produces alcohol + acetate, NO amine)
        acetate_chebis = {
            "CHEBI:30089",  # acetate (from N-acetyl compounds)
            "CHEBI:40480",  # acetate (another form)
        }

        # Products indicating other amide cleavages
        other_amide_products = {
            "CHEBI:16199",  # urea (from guanidino compounds)
            "CHEBI:57845",  # urea (charged form)
            "CHEBI:57972",  # beta-alanine (from pantothenate)
            "CHEBI:16958",  # beta-alanine (another form)
            "CHEBI:57305",  # glycine (from N-acylglycines)
            "CHEBI:15428",  # glycine (another form)
            "CHEBI:29947",  # formate (from formamides)
            "CHEBI:15740",  # formate
        }

        has_nitrogen_product = any(
            p.chebi_id in nitrogen_products
            for p in reaction.right_participants
        )

        has_acetate_product = any(
            p.chebi_id in acetate_chebis
            for p in reaction.right_participants
        )

        has_other_amide_product = any(
            p.chebi_id in other_amide_products
            for p in reaction.right_participants
        )

        # Acetate release only counts as amidase if:
        # 1. There's also an amine/ammonia product (true N-deacetylation), OR
        # 2. The substrate has an explicit amide motif
        has_deacetylation = has_acetate_product and (has_nitrogen_product or has_amide_substrate)

        # Must have cleavage-product evidence that is characteristic of C-N hydrolysis.
        has_cleavage_product_evidence = (
            has_nitrogen_product or has_other_amide_product or has_deacetylation
        )
        if not has_cleavage_product_evidence:
            return ClassificationResult(
                is_member=False,
                explanation="No amide bond cleavage pattern detected"
            )

        # CRITICAL: Exclude ester hydrolysis (acetyl esters release acetate but no amine)
        # If acetate is produced but NO amine and no amide substrate, it's likely ester hydrolysis
        if has_acetate_product and not has_nitrogen_product and not has_amide_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="Acetate released without amine - likely ester hydrolysis not amidase"
            )

        # Apply exclusions before classifying as amidase

        # Exclude peptidases (act on peptide bonds, not simple amides)
        
        has_peptide = False
        
        if has_peptide:
            return ClassificationResult(
                is_member=False,
                explanation="Peptidase - acts on peptide bonds not simple amides"
            )

        # Exclude nucleotide/nucleoside reactions (different C-N bond type)
        
        has_nucleotide = False
        
        if has_nucleotide:
            return ClassificationResult(
                is_member=False,
                explanation="Nucleotide/nucleoside reaction - not simple amide hydrolysis"
            )

        # Exclude complex biosynthetic reactions (ATP involvement)
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        
        if has_atp:
            return ClassificationResult(
                is_member=False,
                explanation="ATP-dependent reaction - not simple amide hydrolysis"
            )

        # Exclude transferase reactions (group transfer, not hydrolysis)
        has_coa = False
        
        # Count non-water substrates
        non_water_substrates = [
            p for p in reaction.left_participants
            if p.chebi_id != CHEBI_H2O]
        
        if has_coa and len(non_water_substrates) > 1:
            return ClassificationResult(
                is_member=False,
                explanation="CoA-dependent transferase - not simple amide hydrolysis"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Amidase: amide bond hydrolysis"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)

    @classmethod
    def _has_amide_moiety(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        mol = participant.get_mol()
        if mol is None or cls.AMIDE_PATTERN is None:
            return False
        return mol.HasSubstructMatch(cls.AMIDE_PATTERN)
