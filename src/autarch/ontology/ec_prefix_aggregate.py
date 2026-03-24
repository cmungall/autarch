"""Helpers for aggregate reaction classifiers.

These utilities implement coarse EC parent classes as aggregates over more
specific curated children already present in the ontology. This keeps the broad
parent definitions parsimonious and transparent: if a parent fires, it is
because one of its more specific children fired.
"""

from __future__ import annotations

from functools import lru_cache
from typing import ClassVar, Type, cast

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass


def _normalize_ec(ec: str) -> str:
    parts = ec.replace("EC:", "").split(".")
    while len(parts) < 4:
        parts.append("-")
    return ".".join(parts[:4])


def _ec_specificity(ec: str) -> int:
    return sum(1 for part in _normalize_ec(ec).split(".") if part != "-")


def _ec_matches_prefix(ec_number: str, ec_prefix: str) -> bool:
    ec_number_parts = _normalize_ec(ec_number).split(".")
    ec_prefix_parts = _normalize_ec(ec_prefix).split(".")
    for number_part, prefix_part in zip(ec_number_parts, ec_prefix_parts):
        if prefix_part == "-":
            return True
        if number_part != prefix_part:
            return False
    return True


def _all_reaction_classes() -> tuple[Type[ReactionClass], ...]:
    import autarch.ontology  # noqa: F401

    def get_all_subclasses(base: type[object]) -> set[type[object]]:
        subclasses: set[type[object]] = set()
        for subclass in base.__subclasses__():
            subclasses.add(subclass)
            subclasses.update(get_all_subclasses(subclass))
        return subclasses

    return tuple(
        cast(Type[ReactionClass], cls)
        for cls in get_all_subclasses(ReactionClass)
        if issubclass(cls, ReactionClass)
        if not getattr(cls, "__abstractmethods__", None)
        if not cls.__dict__.get("EXCLUDE_FROM_DISCOVERY", False)
    )


def _class_ec_prefix(cls: Type[ReactionClass]) -> str | None:
    ec_prefix = getattr(cls, "EC_NUMBER_PREFIX", None)
    if isinstance(ec_prefix, str):
        return ec_prefix
    ec_number = getattr(cls, "EC_NUMBER", None)
    if isinstance(ec_number, str):
        return ec_number
    return None


@lru_cache(maxsize=None)
def aggregate_child_classes(
    ec_prefix: str,
    current_class_name: str,
) -> tuple[Type[ReactionClass], ...]:
    parent_specificity = _ec_specificity(ec_prefix)
    matching_children: list[Type[ReactionClass]] = []
    for cls in _all_reaction_classes():
        if cls.__name__ == current_class_name:
            continue
        child_prefix = _class_ec_prefix(cls)
        if child_prefix is None:
            continue
        if not _ec_matches_prefix(child_prefix, ec_prefix):
            continue
        if _ec_specificity(child_prefix) <= parent_specificity:
            continue
        matching_children.append(cls)

    matching_children.sort(
        key=lambda cls: (-_ec_specificity(_class_ec_prefix(cls) or ""), cls.__name__)
    )
    return tuple(matching_children)


def aggregate_ec_prefix_membership(
    reaction: Reaction,
    current_cls: Type[ReactionClass],
    ec_prefix: str,
    label: str,
) -> ClassificationResult:
    for child_cls in aggregate_child_classes(ec_prefix, current_cls.__name__):
        result = child_cls().check_membership(reaction)
        if result.is_member:
            child_prefix = _class_ec_prefix(child_cls) or "unknown EC"
            return ClassificationResult(
                is_member=True,
                explanation=f"{label}: matched {child_cls.__name__} ({child_prefix}); {result.explanation}",
            )

    return ClassificationResult(
        is_member=False,
        explanation=f"No supported {label.lower()} child classifier matched",
    )


def aggregate_ec_prefix_supports_evaluation(
    reaction: Reaction,
    current_cls: Type[ReactionClass],
    ec_prefix: str,
) -> bool:
    children = aggregate_child_classes(ec_prefix, current_cls.__name__)
    return any(child_cls.supports_evaluation(reaction) for child_cls in children)


def aggregate_explicit_child_membership(
    reaction: Reaction,
    child_classes: tuple[Type[ReactionClass], ...],
    label: str,
) -> ClassificationResult:
    """Aggregate membership over an explicit child-class list."""
    for child_cls in child_classes:
        result = child_cls().check_membership(reaction)
        if result.is_member:
            child_prefix = _class_ec_prefix(child_cls) or "unknown EC"
            return ClassificationResult(
                is_member=True,
                explanation=f"{label}: matched {child_cls.__name__} ({child_prefix}); {result.explanation}",
            )

    return ClassificationResult(
        is_member=False,
        explanation=f"No supported {label.lower()} child classifier matched",
    )


def aggregate_explicit_child_supports_evaluation(
    reaction: Reaction,
    child_classes: tuple[Type[ReactionClass], ...],
) -> bool:
    """Return True if any explicit child can evaluate the reaction."""
    return any(child_cls.supports_evaluation(reaction) for child_cls in child_classes)


class ExplicitEcAggregate(ReactionClass):
    """Base class for explicit EC aggregate wrappers over child classifiers."""

    EXCLUDE_FROM_DISCOVERY: ClassVar[bool] = True
    CHILD_CLASSES: ClassVar[tuple[type[ReactionClass], ...]]

    @classmethod
    def concept_phrase(cls) -> str:
        """Return a human-readable label for aggregate explanations."""
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
            type(self).CHILD_CLASSES,
            type(self).concept_phrase(),
        )
