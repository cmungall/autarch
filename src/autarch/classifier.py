"""Reaction classification convenience methods.

This module provides a high-level interface for classifying reactions using
all available ReactionClass implementations.
"""

from typing import Type

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology import ReactionClass


class ReactionClassifier:
    """Main classifier that can use all registered reaction classes.

    >>> from autarch.classifier import ReactionClassifier
    >>> from autarch.ontology import Reaction, Participant
    >>> classifier = ReactionClassifier()

    >>> # Test ester hydrolysis: ethyl acetate + water → acetic acid + ethanol
    >>> reaction = Reaction(
    ...     left_participants=[
    ...         Participant(smiles="CC(=O)OCC"),  # ethyl acetate
    ...         Participant(smiles="O")  # water
    ...     ],
    ...     right_participants=[
    ...         Participant(smiles="CC(=O)O"),  # acetic acid
    ...         Participant(smiles="CCO")  # ethanol
    ...     ]
    ... )
    >>> results = classifier.classify(reaction)
    >>> 'Hydrolase' in results
    True
    >>> results['Hydrolase'].is_member
    True

    >>> # Test with specific class
    >>> from autarch.ontology import Hydrolase
    >>> result = classifier.classify_with_class(reaction, Hydrolase)
    >>> result.is_member
    True

    >>> # Test non-hydrolysis reaction: combustion
    >>> combustion = Reaction(
    ...     left_participants=[
    ...         Participant(smiles="CC"),  # ethane
    ...         Participant(smiles="O=O", count=7)  # oxygen
    ...     ],
    ...     right_participants=[
    ...         Participant(smiles="O=C=O", count=2),  # CO2
    ...         Participant(smiles="O", count=3)  # water
    ...     ]
    ... )
    >>> results = classifier.classify(combustion)
    >>> results.get('Hydrolase').is_member if 'Hydrolase' in results else False
    False
    """

    def __init__(self):
        """Initialize the classifier and discover all available reaction classes."""
        self.reaction_classes = self._discover_reaction_classes()

    def _discover_reaction_classes(self) -> dict[str, Type[ReactionClass]]:
        """Discover all subclasses of ReactionClass.

        Returns:
            Dictionary mapping class names to class types
        """

        def get_all_subclasses(cls):
            """Recursively get all subclasses of a class."""
            subclasses = set()
            for subclass in cls.__subclasses__():
                subclasses.add(subclass)
                subclasses.update(get_all_subclasses(subclass))
            return subclasses

        # Get all concrete subclasses (not abstract)
        all_subclasses = get_all_subclasses(ReactionClass)
        concrete_classes = {
            cls.__name__: cls
            for cls in all_subclasses
            if not getattr(cls, "__abstractmethods__", None)
            if not cls.__dict__.get("EXCLUDE_FROM_DISCOVERY", False)
        }

        return concrete_classes

    def classify(self, reaction: Reaction) -> dict[str, ClassificationResult]:
        """Classify a reaction against all known reaction classes.

        Args:
            reaction: The reaction to classify

        Returns:
            Dictionary mapping class names to classification results
        """
        results = {}
        for class_name, reaction_class in self.reaction_classes.items():
            instance = reaction_class()
            result = instance.check_membership(reaction)
            results[class_name] = result

        return results

    def classify_with_class(
        self, reaction: Reaction, reaction_class: Type[ReactionClass]
    ) -> ClassificationResult:
        """Classify a reaction using a specific reaction class.

        Args:
            reaction: The reaction to classify
            reaction_class: The specific ReactionClass to use

        Returns:
            Classification result
        """
        instance = reaction_class()
        return instance.check_membership(reaction)

    def get_matching_classes(self, reaction: Reaction) -> list[str]:
        """Get names of all reaction classes that match the given reaction.

        Args:
            reaction: The reaction to classify

        Returns:
            List of class names where the reaction is a member
        """
        results = self.classify(reaction)
        return [
            class_name for class_name, result in results.items() if result.is_member
        ]


def classify_reaction(reaction: Reaction) -> dict[str, ClassificationResult]:
    """Convenience function to classify a reaction.

    >>> from autarch.classifier import classify_reaction
    >>> from autarch.ontology import Reaction, Participant

    >>> # Simple ester hydrolysis
    >>> reaction = Reaction(
    ...     left_participants=[
    ...         Participant(smiles="CC(=O)OC"),  # methyl acetate
    ...         Participant(smiles="O")  # water
    ...     ],
    ...     right_participants=[
    ...         Participant(smiles="CC(=O)O"),  # acetic acid
    ...         Participant(smiles="CO")  # methanol
    ...     ]
    ... )
    >>> results = classify_reaction(reaction)
    >>> 'Hydrolase' in results
    True
    """
    classifier = ReactionClassifier()
    return classifier.classify(reaction)


def get_reaction_classes() -> list[str]:
    """Get names of all available reaction classes.

    >>> from autarch.classifier import get_reaction_classes
    >>> classes = get_reaction_classes()
    >>> 'Hydrolase' in classes
    True
    """
    classifier = ReactionClassifier()
    return list(classifier.reaction_classes.keys())
