"""Helpers for declarative GO aggregate classifiers."""

from __future__ import annotations

from typing import ClassVar

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.reaction import ReactionClass


class ExplicitGoAggregate(ReactionClass):
    """Base class for GO-backed wrappers over explicit child classifiers."""

    EXCLUDE_FROM_DISCOVERY: ClassVar[bool] = True
    CHILD_CLASSES: ClassVar[tuple[type[ReactionClass], ...]]

    @classmethod
    def concept_phrase(cls) -> str:
        """Return a human-readable label for aggregate explanations.

        Prefer an explicitly declared ``CONCEPT_PHRASE`` when present, but
        default to the class docstring so small union-of wrappers do not need
        extra boilerplate.
        """
        explicit = getattr(cls, "CONCEPT_PHRASE", None)
        if isinstance(explicit, str) and explicit:
            return explicit

        docstring = (cls.__doc__ or "").strip()
        if docstring:
            return docstring.splitlines()[0].strip().rstrip(".")

        return cls.__name__

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            type(self).concept_phrase(),
        )
