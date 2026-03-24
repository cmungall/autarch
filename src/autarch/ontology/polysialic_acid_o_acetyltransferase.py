"""Polysialic Acid O-Acetyltransferase reaction classification.

Polysialic acid O-acetyltransferase (EC 2.3.1.136) catalyzes the O-acetylation
of polysialic acid chains:

    [N-acetyl-α-D-neuraminosyl-(2->8)](n) + n acetyl-CoA →
    [O-acetyl-N-acetyl-α-D-neuraminosyl-(2->8)](n) + n CoA

This modifies sialic acid polymers by adding O-acetyl groups.
"""

from autarch.datamodel import ClassificationResult, Reaction, PolymerType
from autarch.ontology.reaction import IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE, ReactionClass


class PolysialicAcidOAcetyltransferase(ReactionClass):
    """Classifier for polysialic-acid O-acetyltransferase reactions.

    Detects O-acetylation of polysialic acid polymers:
    - Neuraminosyl/sialic acid polymer on both sides
    - Acetyl-CoA as acetyl donor
    - CoA as leaving group
    - Addition of O-acetyl group to the polymer
    """

    GO_ID = "GO:0050208"  # polysialic-acid O-acetyltransferase activity
    EC_NUMBERS = ["2.3.1.136"]
    EC_NUMBER_PREFIX = None
    EVALUATION_EVIDENCE = IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE
    ACETYL_COA_CHEBI = "CHEBI:57288"
    COA_CHEBI = "CHEBI:57287"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is catalyzed by polysialic acid O-acetyltransferase.

        Detects:
        - Sialic acid polymer on both sides
        - Acetyl-CoA/CoA involvement
        - O-acetylation pattern
        """
        # Check for sialic acid polymer
        left_polymer = self._find_sialic_acid_polymer(reaction.left_participants)
        right_polymer = self._find_sialic_acid_polymer(reaction.right_participants)

        if not left_polymer or not right_polymer:
            return ClassificationResult(
                is_member=False,
                explanation="No sialic acid polymer found on both sides"
            )

        left_acetyl_coa = self._has_chebi(reaction.left_participants, self.ACETYL_COA_CHEBI)
        right_acetyl_coa = self._has_chebi(reaction.right_participants, self.ACETYL_COA_CHEBI)
        left_coa = self._has_chebi(reaction.left_participants, self.COA_CHEBI)
        right_coa = self._has_chebi(reaction.right_participants, self.COA_CHEBI)

        if not ((left_acetyl_coa or right_acetyl_coa) and (left_coa or right_coa)):
            return ClassificationResult(
                is_member=False,
                explanation="No acetyl-CoA/CoA pair found"
            )

        if left_polymer.polymer_index and right_polymer.polymer_index:
            if left_polymer.polymer_index != right_polymer.polymer_index:
                return ClassificationResult(
                    is_member=False,
                    explanation="Polymer stoichiometry/index does not align across the transfer",
                )

        if left_acetyl_coa and right_coa and not right_acetyl_coa and not left_coa:
            explanation = "Polysialic acid O-acetyltransferase: acetyl-CoA-dependent O-acetylation of a sialic acid polymer"
        elif left_coa and right_acetyl_coa and not left_acetyl_coa and not right_coa:
            explanation = "Polysialic acid O-acetyltransferase (reverse): deacetylation of an O-acetylated sialic acid polymer"
        else:
            return ClassificationResult(
                is_member=False,
                explanation="Acetyl-CoA/CoA direction does not support polysialic acid O-acetyl transfer"
            )

        return ClassificationResult(
            is_member=True,
            explanation=explanation
        )

    def _find_sialic_acid_polymer(self, participants):
        """Find sialic acid/neuraminosyl polymer participant."""
        for p in participants:
            if p.polymer_type == PolymerType.SIALIC_ACID_POLYMER:
                return p
        return None

    @staticmethod
    def _has_chebi(participants, chebi_id: str) -> bool:
        """Check if a participant list contains the given ChEBI identifier."""
        return any(participant.chebi_id == chebi_id for participant in participants)
