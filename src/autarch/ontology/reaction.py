from abc import ABC, abstractmethod
from enum import Enum
from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction


class EvaluationEvidence(str, Enum):
    """Types of evidence sufficient to evaluate a reaction for one classifier."""

    COMPLETE_SMILES = "complete_smiles"
    POLYMER_METADATA = "polymer_metadata"
    PARTICIPANT_IDENTIFIERS = "participant_identifiers"


IDENTIFIER_FRIENDLY_EVIDENCE = frozenset(
    {
        EvaluationEvidence.COMPLETE_SMILES,
        EvaluationEvidence.PARTICIPANT_IDENTIFIERS,
    }
)
"""Evidence modes for classifiers that can run from participant IDs alone."""


IDENTIFIER_AND_POLYMER_FRIENDLY_EVIDENCE = frozenset(
    {
        EvaluationEvidence.COMPLETE_SMILES,
        EvaluationEvidence.POLYMER_METADATA,
        EvaluationEvidence.PARTICIPANT_IDENTIFIERS,
    }
)
"""Evidence modes for classifiers that can run from IDs or polymer metadata."""


class ReactionClass(ABC):
    """Abstract base class for reaction classifiers.

    All subclasses must define GO_ID with the appropriate GO term.
    """

    GO_ID: ClassVar[Optional[str]] = None  # Subclasses must override with GO term ID
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = None  # Optional EC number prefix
    EC_BROAD_XREFS: ClassVar[list[str]] = []  # Optional non-exact EC cross-references
    EVALUATION_EVIDENCE: ClassVar[frozenset[EvaluationEvidence]] = frozenset(
        {
            EvaluationEvidence.COMPLETE_SMILES,
            EvaluationEvidence.POLYMER_METADATA,
        }
    )

    def check_membership(self, reaction: Reaction) -> ClassificationResult:
        """Check if a given reaction is a member of the reaction class."""
        # Ensure all molecules are parsed
        for stoi in reaction.all_participants():
            stoi.get_mol()
        return self.check_membership_impl(reaction)

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        """Return True if the reaction carries enough evidence for this classifier."""
        participants = reaction.all_participants()

        if (
            EvaluationEvidence.COMPLETE_SMILES in cls.EVALUATION_EVIDENCE
            and all(participant.smiles is not None for participant in participants)
        ):
            return True

        if (
            EvaluationEvidence.POLYMER_METADATA in cls.EVALUATION_EVIDENCE
            and any(
                participant.polymer_index is not None
                or participant.polymer_type is not None
                for participant in participants
            )
        ):
            return True

        if (
            EvaluationEvidence.PARTICIPANT_IDENTIFIERS in cls.EVALUATION_EVIDENCE
            and all(
                participant.chebi_id is not None
                or participant.name is not None
                or participant.polymer_index is not None
                for participant in participants
            )
        ):
            return True

        return False

    @abstractmethod
    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if a given RDKit Mol object is a member of the reaction class."""
        pass
