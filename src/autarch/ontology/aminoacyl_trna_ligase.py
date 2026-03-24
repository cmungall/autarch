"""tRNA Synthetase reaction classification using pattern DSL.

tRNA Synthetases are aminoacyl-tRNA synthetases that attach amino acids to their
corresponding tRNA molecules. This is a challenging class because it involves:
1. Polymer substrates (tRNA)
2. Very specific amino acid recognition
3. ATP-dependent aminoacyl activation

Key biochemical signature: amino acid + tRNA + ATP → aminoacyl-tRNA + AMP + PPi

This is one of the most ancient and conserved enzymatic reactions, essential
for protein synthesis.

History
-------

## 2025-12-22

Fixed false positives by:
1. Requiring "generic polymer" participant (tRNA indicator in RHEA)
2. Excluding CoA-containing reactions (CoA ligases, not tRNA synthetases)
3. Excluding asparagine/glutamine synthetases (produce amino acids, not aminoacyl-tRNA)
Previously matched ANY ATP→AMP+PPi reaction with 3→3/3→4 stoichiometry.
"""

from autarch.datamodel import ClassificationResult, PolymerType, Reaction
from autarch.ontology.ligase import Ligase
from autarch.pattern_dsl import var, optional
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    amp,
    atp,
    diphosphate,
    h_plus,
    p,
)


class AminoacylTRNALigase(Ligase):
    """aminoacyl tRNA ligase

    Examples:
    - Tyrosyl-tRNA synthetase: tRNA(Tyr) + Tyr + ATP → Tyr-tRNA(Tyr) + AMP + PPi
    - Valyl-tRNA synthetase: tRNA(Val) + Val + ATP → Val-tRNA(Val) + AMP + PPi
    - Leucyl-tRNA synthetase: tRNA(Leu) + Leu + ATP → Leu-tRNA(Leu) + AMP + PPi

    This class tests our ability to handle:
    - Polymer substrates (tRNA molecules)
    - Naming patterns (aminoacyl-tRNA formation)
    - ATP → AMP + PPi energetics (different from kinases)
    """

    # NOTE: GO:0004812 (aminoacyl-tRNA ligase activity) not present in evaluation dataset
    # Using EC-only evaluation since EC 6.1.1.* is the true equivalence
    GO_ID = "GO:0004812"  # aminoacyl-tRNA ligase activity
    EC_NUMBER_PREFIX = "6.1.1.-"  # Ligases forming aminoacyl-tRNA

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Classic aminoacyl-tRNA synthetase: tRNA + amino_acid + ATP → aminoacyl-tRNA + AMP + PPi
        var("tRNA") + var("amino_acid") + p(atp) 
        >> var("aminoacyl_tRNA") + p(amp) + p(diphosphate) + optional(h_plus),
        # Without explicit H+
        var("tRNA") + var("amino_acid") + p(atp) 
        >> var("aminoacyl_tRNA") + p(amp) + p(diphosphate),
        # Alternative order (some reactions have different participant order)
        var("amino_acid") + var("tRNA") + p(atp)
        >> var("aminoacyl_tRNA") + p(amp) + p(diphosphate),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a tRNA synthetase using SIMPLE, DECLARATIVE approach.

        BIOCHEMICAL SIGNATURE:
        1. Must have ATP as energy source
        2. Must produce AMP + PPi (not ADP + Pi like kinases)
        3. Must involve tRNA (polymer substrate - "generic polymer" in RHEA)
        4. Must NOT be CoA ligase, asparagine synthetase, or other non-tRNA ligases
        5. Should be exactly 3→3 or 3→4 reaction (ATP activation pattern)
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must have ATP as energy source
        has_atp = any(
            p.chebi_id == CHEBI_ATP
            for p in reaction.left_participants
        )

        if not has_atp:
            return ClassificationResult(
                is_member=False,
                explanation="No ATP - tRNA synthetases require ATP for aminoacyl activation"
            )

        # Must produce AMP (not ADP) - this distinguishes from kinases
        has_amp = any(
            p.chebi_id == CHEBI_AMP
            for p in reaction.right_participants
        )

        if not has_amp:
            return ClassificationResult(
                is_member=False,
                explanation="No AMP product - tRNA synthetases use ATP → AMP + PPi energetics"
            )

        # Must produce diphosphate (PPi)
        has_diphosphate = any(
            p.chebi_id == CHEBI_DIPHOSPHATE
            for p in reaction.right_participants
        )

        if not has_diphosphate:
            return ClassificationResult(
                is_member=False,
                explanation="No diphosphate - tRNA synthetases produce AMP + PPi"
            )

        # CRITICAL: Must involve tRNA - now using polymer_type from ETL
        # tRNA synthetases have tRNA on BOTH sides (tRNA → aminoacyl-tRNA)
        has_trna_substrate = any(
            part.polymer_type == PolymerType.TRNA
            for part in reaction.left_participants
        )

        has_trna_product = any(
            part.polymer_type == PolymerType.TRNA
            for part in reaction.right_participants
        )

        if not (has_trna_substrate and has_trna_product):
            return ClassificationResult(
                is_member=False,
                explanation="No tRNA on both sides - not aminoacyl-tRNA synthesis"
            )

        # Exclude CoA ligases (have CoA as substrate or CoA in product name)
        CHEBI_COA = "CHEBI:57287"  # coenzyme A(4-)
        has_coa = any(
            p.chebi_id == CHEBI_COA
            for p in reaction.left_participants + reaction.right_participants
        ) or "coa" in label_lower

        if has_coa:
            return ClassificationResult(
                is_member=False,
                explanation="CoA present - CoA ligase not tRNA synthetase"
            )

        # Exclude asparagine/glutamine synthetases (produce amino acids, not aminoacyl-tRNA)
        # These produce asparagine/glutamine + AMP + PPi, not charged tRNA
        AMINO_ACID_PRODUCTS = {"asparagine", "glutamine", "pantothenate"}
        # Check if any of these appear in product_str outside the polymer context
        has_amino_acid_product = any(aa in product_str for aa in AMINO_ACID_PRODUCTS) and "polymer" not in product_str

        if has_amino_acid_product:
            return ClassificationResult(
                is_member=False,
                explanation="Produces free amino acid - amino acid synthetase not tRNA synthetase"
            )

        # Should be 3→3, 3→4, or 3→5 reaction
        # (tRNA + amino acid + ATP → charged tRNA + AMP + PPi [+ H+] [+ duplicates from ETL])
        # Note: ETL may add duplicate participants during label matching
        reactant_count = len(reaction.left_participants)
        product_count = len(reaction.right_participants)

        if reactant_count != 3 or product_count not in [3, 4, 5]:
            return ClassificationResult(
                is_member=False,
                explanation=f"Wrong stoichiometry ({reactant_count}→{product_count}) - expected 3→3/3→4/3→5"
            )

        return ClassificationResult(
            is_member=True,
            explanation="Aminoacyl tRNA ligase: aminoacyl-tRNA formation (polymer + AA + ATP -> charged polymer + AMP + PPi)"
        )
