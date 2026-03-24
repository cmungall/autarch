"""Ligase forming phosphoric ester bonds classification.

EC 6.5 ligases catalyze the formation of phosphodiester bonds,
primarily DNA and RNA ligases that seal nicks in nucleic acid strands.
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ligase import Ligase


class LigaseFormingPhosphoricEsterBonds(Ligase):
    """ligase activity, forming phosphoric ester bonds

    Catalysis of the joining of two molecules via a phosphoric ester bond,
    with the concomitant hydrolysis of ATP or a similar triphosphate.

    Examples:
    - DNA ligase: seals single-strand nicks in double-stranded DNA
    - RNA ligase: joins RNA fragments
    - Polynucleotide ligase

    These enzymes form phosphodiester bonds in nucleic acid backbones.
    """

    GO_ID = "GO:0016886"  # ligase activity, forming phosphoric ester bonds
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "6.5.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a phosphoric ester bond-forming ligase.

        Strategy:
        1. Must be a ligase (ATP-dependent bond formation)
        2. Look for DNA/RNA/polynucleotide substrates by label
        """
        # First check if it's a ligase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a ligase: {parent_result.explanation}",
            )

        # Label-based detection for nucleic acid ligases
        all_names = [
            p.name.lower()
            for p in reaction.left_participants + reaction.right_participants
            if p.name
        ]
        nucleic_acid_terms = {
            "dna ligase", "rna ligase", "polynucleotide",
            "dna nick", "phosphodiester",
        }
        if any(term in name for name in all_names for term in nucleic_acid_terms):
            return ClassificationResult(
                is_member=True,
                explanation="Phosphoric ester ligase: DNA/RNA ligase detected by label",
            )

        # Check for DNA/RNA as participants (broad label check)
        dna_rna_terms = {"dna", "rna", "nucleotide"}
        has_nucleic_acid = any(
            any(term in name for term in dna_rna_terms)
            for name in all_names
        )
        if has_nucleic_acid:
            return ClassificationResult(
                is_member=True,
                explanation="Phosphoric ester ligase: nucleic acid substrate with ATP-dependent ligation",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No nucleic acid substrates - not a phosphoric ester bond-forming ligase",
        )
