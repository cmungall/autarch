"""Oxidoreductase acting on CH or CH2 groups reaction classification.

EC 1.17: Oxidoreductases acting on CH or CH2 groups.
These enzymes oxidize C-H bonds in contexts that are NOT alcohol (C-OH),
amine (C-NH2/C-NH), or aldehyde (C=O) groups. Examples include xanthine
oxidase, ribonucleoside-diphosphate reductase, and 4-hydroxyphenylpyruvate
dioxygenase.

Key features:
- Disulfide/thioredoxin acceptors (ribonucleoside-diphosphate reductase)
- O2 acceptors with specialized substrates (xanthine oxidase)
- Oxidation at sp3 carbon in non-standard contexts
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase

# Thioredoxin (disulfide/dithiol) - key acceptor for EC 1.17.4
CHEBI_THIOREDOXIN_DISULFIDE = "CHEBI:50058"  # oxidized thioredoxin
CHEBI_THIOREDOXIN_DITHIOL = "CHEBI:29950"  # reduced thioredoxin


class OxidoreductaseActingOnCHOrCH2Groups(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction in which a CH or CH2 group acts as hydrogen or electron donor and reduces an acceptor."""

    GO_ID = "GO:0016725"  # oxidoreductase activity, acting on CH or CH2 groups
    EC_NUMBER_PREFIX = "1.17.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase acting on CH or CH2 groups.

        Strategy:
        1. Must be an oxidoreductase (parent class check)
        2. Look for thioredoxin/disulfide acceptors (EC 1.17.4 - reductases)
        3. Look for label indicators: xanthine, ribonucleoside-diphosphate reductase
        """
        # First check parent
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not an oxidoreductase: {parent_result.explanation}",
            )

        label_lower = reaction.label.lower() if reaction.label else ""

        # Check for thioredoxin as acceptor (EC 1.17.4 - ribonucleoside-diphosphate reductase)
        has_thioredoxin = any(
            p.chebi_id in {CHEBI_THIOREDOXIN_DISULFIDE, CHEBI_THIOREDOXIN_DITHIOL}
            for p in reaction.left_participants + reaction.right_participants
        ) or "thioredoxin" in label_lower

        # Label-based detection for well-known EC 1.17 reactions
        ch_oxidation_labels = [
            "xanthine",
            "ribonucleoside-diphosphate reductase",
            "ribonucleoside diphosphate reductase",
            "4-hydroxyphenylpyruvate",
        ]
        has_ch_label = any(pat in label_lower for pat in ch_oxidation_labels)

        if has_thioredoxin:
            return ClassificationResult(
                is_member=True,
                explanation="CH/CH2 oxidoreductase: thioredoxin-dependent reductase",
            )

        if has_ch_label:
            return ClassificationResult(
                is_member=True,
                explanation=f"CH/CH2 oxidoreductase: label match in '{label_lower}'",
            )

        return ClassificationResult(
            is_member=False,
            explanation="No CH/CH2 oxidoreductase pattern detected",
        )
