"""amino acid transmembrane transporter activity."""

from __future__ import annotations

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.transport_utils import reactive_transported_pairs


class AminoAcidTransmembraneTransporterActivity(ReactionClass):
    """amino acid transmembrane transporter activity."""

    GO_ID = "GO:0015171"
    ALPHA_AMINO_ACID_PATTERN = Chem.MolFromSmarts(
        "[NX3,NX4+;!$([N+](C)(C)C)][CX4][CX3](=[OX1])[O-,$([OX2H1])]"
    )

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for amino-acid transport across a membrane."""
        if not reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Not a transport reaction",
            )

        transported = reactive_transported_pairs(reaction)
        if not transported:
            return ClassificationResult(
                is_member=False,
                explanation="No transported substrate pair beyond transport spectators",
            )

        if not all(self._is_amino_acid_like(left) for left, _ in transported):
            return ClassificationResult(
                is_member=False,
                explanation="Transported substrate is not consistently amino-acid-like",
            )

        transport_info = [
            f"{left.chebi_id or 'unknown'}: {left.location or 'unknown'}→{right.location or 'unknown'}"
            for left, right in transported
        ]
        return ClassificationResult(
            is_member=True,
            explanation=f"Amino acid transmembrane transporter: {', '.join(transport_info)}",
        )

    @classmethod
    def _is_amino_acid_like(cls, participant: Participant) -> bool:
        mol = participant.get_mol()
        if mol is None or cls.ALPHA_AMINO_ACID_PATTERN is None:
            return False
        return mol.HasSubstructMatch(cls.ALPHA_AMINO_ACID_PATTERN)
