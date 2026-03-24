"""
Aminopeptidase classification - EC 3.4.11.

Aminopeptidases are hydrolases that cleave amino acids from the N-terminal end of peptides:

N-terminal-amino-acid-peptide + H2O → amino acid + peptide

Examples:
- Leucyl aminopeptidase (EC 3.4.11.1)
- Lysyl aminopeptidase (EC 3.4.11.15)
- Arginyl aminopeptidase (EC 3.4.11.27)
- Tryptophanyl aminopeptidase (EC 3.4.11.17)
"""

from rdkit import Chem

from autarch.ontology.reaction import ReactionClass
from autarch.datamodel import Reaction, ClassificationResult

from autarch.molecules import CHEBI_H2O


# SMARTS pattern for peptide bond: C(=O)N-C
PEPTIDE_BOND_SMARTS = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")

# L-amino acid backbone pattern
AMINO_ACID_SMARTS = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")


class Aminopeptidase(ReactionClass):
    """aminopeptidase

    Detects cleavage of N-terminal amino acids from peptides.
    """

    GO_ID = "GO:0004177"  # aminopeptidase activity
    EC_NUMBER_PREFIX = "3.4.11.-"  # EC prefix for aminopeptidases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """
        Check if reaction represents aminopeptidase activity.

        Strategy:
        1. Look for water as reactant
        2. Detect peptide substrate (multiple peptide bonds)
        3. Look for amino acid products
        4. Verify shorter peptide product (fewer peptide bonds)
        """
        # Use RHEA label for pattern matching
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Check for water reactant
        has_water = any(
            p.chebi_id == CHEBI_H2O or (p.smiles == "O")
            for p in reaction.left_participants
        )

        if not has_water:
            return ClassificationResult(
                is_member=False, explanation="No water reactant found"
            )

        # Count peptide bonds in reactants
        # Note: A single "peptide bond" SMARTS match isn't enough because it also
        # matches cyclic amides (e.g., 5-oxoproline). We require at least 2 matches
        # to indicate a true peptide chain, OR explicit naming.
        reactant_peptide_bonds = 0
        has_peptide_substrate = False

        for p in reaction.left_participants:
            mol = p.get_mol()
            if mol is not None and PEPTIDE_BOND_SMARTS is not None:
                matches = mol.GetSubstructMatches(PEPTIDE_BOND_SMARTS)
                if matches:
                    reactant_peptide_bonds += len(matches)
                    # Require at least 2 peptide bonds for structural detection
                    # (dipeptide = 1 bond, but we can't distinguish from cyclic amides)
                    if len(matches) >= 2:
                        has_peptide_substrate = True

        if not has_peptide_substrate:
            # Check for peptide-related patterns in label (more reliable for dipeptides)
            # Use specific patterns: "peptide", "polypeptide", "protein", or
            # amino acid dipeptide naming like "L-alanyl-L-glycine" (amino-yl-amino)
            # Explicit peptide/protein keywords
            if any(kw in substrate_str for kw in ["peptide", "polypeptide", "protein"]):
                has_peptide_substrate = True
                reactant_peptide_bonds = max(reactant_peptide_bonds, 1)
            # Dipeptide pattern: check for amino-yl-amino naming
            # e.g., "L-leucyl-L-proline" but NOT "N-acetyl-D-glucosaminyl"
            elif "-yl-" in substrate_str and "l-" in substrate_str:
                # Exclude glycosyl/acetyl compounds
                if not any(exc in substrate_str for exc in ["glucosaminyl", "acetyl", "glycosyl", "mannosyl", "galactosyl"]):
                    has_peptide_substrate = True
                    reactant_peptide_bonds = max(reactant_peptide_bonds, 1)

        if not has_peptide_substrate:
            return ClassificationResult(
                is_member=False, explanation="No peptide substrate found"
            )

        # Check for amino acid products
        # Known amino acids that can be cleaved by aminopeptidases
        AMINO_ACID_CHEBI = {
            "CHEBI:16449",  # alanine
            "CHEBI:29016",  # L-alanine
            "CHEBI:16467",  # arginine
            "CHEBI:32682",  # L-arginine
            "CHEBI:17196",  # asparagine
            "CHEBI:17053",  # aspartic acid
            "CHEBI:15356",  # cysteine
            "CHEBI:16015",  # glutamic acid
            "CHEBI:16017",  # glutamine
            "CHEBI:15428",  # glycine
            "CHEBI:15971",  # histidine
            "CHEBI:16847",  # isoleucine
            "CHEBI:15603",  # leucine
            "CHEBI:25017",  # L-leucine
            "CHEBI:16643",  # lysine
            "CHEBI:16811",  # methionine
            "CHEBI:17295",  # phenylalanine
            "CHEBI:17203",  # proline
            "CHEBI:17115",  # serine
            "CHEBI:16857",  # threonine
            "CHEBI:16828",  # tryptophan
            "CHEBI:17895",  # tyrosine
            "CHEBI:16414",  # valine
            "CHEBI:33704",  # alpha-amino acid
        }

        has_amino_acid_product = False
        for p in reaction.right_participants:
            if p.chebi_id in AMINO_ACID_CHEBI:
                has_amino_acid_product = True
                break
            # Check by structure
            mol = p.get_mol()
            if mol is not None and AMINO_ACID_SMARTS is not None:
                if mol.HasSubstructMatch(AMINO_ACID_SMARTS):
                    # Verify it's a single amino acid (no peptide bonds)
                    peptide_matches = mol.GetSubstructMatches(PEPTIDE_BOND_SMARTS) if PEPTIDE_BOND_SMARTS else []
                    if len(peptide_matches) == 0:
                        has_amino_acid_product = True
                        break

        # Count peptide bonds in products
        product_peptide_bonds = 0
        has_shorter_peptide = False

        for p in reaction.right_participants:
            mol = p.get_mol()
            if mol is not None and PEPTIDE_BOND_SMARTS is not None:
                matches = mol.GetSubstructMatches(PEPTIDE_BOND_SMARTS)
                if matches:
                    product_peptide_bonds += len(matches)
                    has_shorter_peptide = True

        # Exclude glycosyl hydrolases by checking product label
        # Aminopeptidase should not produce sugar/glycan products
        sugar_keywords = ["glucosamin", "galactos", "mannos", "fructos", "xylos", "ribos"]
        has_sugar_product = any(kw in product_str for kw in sugar_keywords)

        if has_sugar_product:
            return ClassificationResult(
                is_member=False,
                explanation="Products include sugars - likely glycosyl hydrolase, not aminopeptidase",
            )

        # Aminopeptidase: peptide substrate + amino acid product + shorter peptide
        if has_peptide_substrate and has_amino_acid_product:
            # For true aminopeptidase:
            # - Substrate has peptide bonds
            # - Product has amino acid + remaining peptide (fewer bonds)
            if product_peptide_bonds < reactant_peptide_bonds and product_peptide_bonds >= 0:
                return ClassificationResult(
                    is_member=True,
                    explanation="Aminopeptidase: N-terminal amino acid cleavage from peptide",
                )
            # Dipeptide case: substrate has 1 peptide bond, products are 2 amino acids
            if reactant_peptide_bonds == 1 and product_peptide_bonds == 0 and has_amino_acid_product:
                return ClassificationResult(
                    is_member=True,
                    explanation="Aminopeptidase: dipeptide → two amino acids",
                )

        # Fallback: peptide → shorter peptide + amino acid
        if has_shorter_peptide and has_amino_acid_product and product_peptide_bonds >= 1:
            return ClassificationResult(
                is_member=True,
                explanation="Aminopeptidase: peptide → amino acid + shorter peptide",
            )

        return ClassificationResult(
            is_member=False, explanation="No clear aminopeptidase pattern detected"
        )
