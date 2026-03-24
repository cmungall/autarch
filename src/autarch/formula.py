"""Chemical reaction formula DSL using operator overloading.

This module provides a natural syntax for writing chemical reactions using
Python's operator overloading capabilities. Since we can't overload the = operator,
we use >> for left-to-right, << for right-to-left, and | for bidirectional reactions.

Basic Usage:
    >>> from autarch.formula import ATP, ADP, Pi, H2O
    >>>
    >>> # Simple ATP hydrolysis reaction
    >>> rxn = ATP + H2O >> ADP + Pi
    >>> print(rxn)
    ATP + H2O → ADP + Pi
    >>>
    >>> # Build to a Reaction object for classification
    >>> reaction = rxn.build()
    >>> len(reaction.left_participants)
    2
    >>> len(reaction.right_participants)
    2
    >>> reaction.left_participants[0].chebi_id
    'CHEBI:15422'

Stoichiometry:
    >>> from autarch.formula import glucose, O2, CO2, H2O
    >>>
    >>> # Cellular respiration with stoichiometry
    >>> respiration = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)
    >>> print(respiration)
    glucose + 6 O2 → 6 CO2 + 6 H2O
    >>>
    >>> # Alternative: left multiplication
    >>> water_formation = (2 * H2O) >> (H * 2) + O2
    >>> reaction = water_formation.build()
    >>> reaction.left_participants[0].count
    2

Compartments/Locations:
    >>> # ATP/ADP translocase with compartment information
    >>> transport = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
    >>> print(transport)
    ATP + ADP → ATP + ADP
    >>>
    >>> # Check it's recognized as transport
    >>> reaction = transport.build()
    >>> reaction.is_transport_reaction()
    True
    >>>
    >>> # Get transport details
    >>> transported = reaction.get_transported_molecules()
    >>> len(transported)
    2
    >>> transported[0][1], transported[0][2]  # from, to locations
    ('out', 'in')

Reaction Directions:
    >>> from autarch.formula import NAD, NADH, H
    >>>
    >>> # Bidirectional reaction using |
    >>> redox = NAD + H | NADH
    >>> redox.direction
    'bidirectional'
    >>> print(redox)
    NAD+ + H+ ⇌ NADH
    >>>
    >>> # Right-to-left using <<
    >>> reverse = ADP + Pi << ATP + H2O
    >>> reverse.direction
    'right-to-left'
    >>> print(reverse)
    ATP + H2O ← ADP + Pi

Custom Molecules:
    >>> from autarch.formula import molecule
    >>>
    >>> # Create custom molecules
    >>> substrate = molecule("CHEBI:99999", name="substrate")
    >>> product = molecule("CHEBI:88888", name="product")
    >>> enzyme_rxn = substrate + ATP >> product + ADP
    >>> print(enzyme_rxn)
    substrate + ATP → product + ADP
    >>>
    >>> # With SMILES instead of CHEBI
    >>> ethanol = molecule(smiles="CCO", name="ethanol")
    >>> acetaldehyde = molecule(smiles="CC=O", name="acetaldehyde")
    >>> oxidation = ethanol + NAD >> acetaldehyde + NADH + H
    >>> print(oxidation)
    ethanol + NAD+ → acetaldehyde + NADH + H+
"""

from __future__ import annotations
from typing import Optional, Union, List
from dataclasses import dataclass, field
from autarch.datamodel import Participant, Reaction


@dataclass
class ParticipantBuilder:
    """Builder for creating Participant objects with operator overloading for reaction formulas.

    The ParticipantBuilder class supports operator overloading for building reactions
    with natural syntax. You can specify stoichiometry, location, and combine
    participants using standard Python operators.

    Examples:
        >>> from autarch.formula import ParticipantBuilder
        >>>
        >>> # Create a participant builder with CHEBI ID
        >>> atp = ParticipantBuilder(chebi_id="CHEBI:15422", name="ATP")
        >>> atp.chebi_id
        'CHEBI:15422'
        >>> atp.name
        'ATP'
        >>> atp.count
        1
        >>>
        >>> # Set stoichiometry with multiplication
        >>> atp2 = atp * 2
        >>> atp2.count
        2
        >>> atp2.name
        'ATP'
        >>>
        >>> # Left multiplication also works
        >>> atp3 = 3 * atp
        >>> atp3.count
        3
        >>>
        >>> # Set location/compartment with subscript
        >>> atp_mito = atp["mitochondria"]
        >>> atp_mito.location
        'mitochondria'
        >>> atp_mito.chebi_id
        'CHEBI:15422'
        >>>
        >>> # Combine stoichiometry and location
        >>> atp_complex = (atp * 2)["cytoplasm"]
        >>> atp_complex.count
        2
        >>> atp_complex.location
        'cytoplasm'
        >>>
        >>> # Add molecules to create ParticipantList
        >>> water = ParticipantBuilder(chebi_id="CHEBI:15377", name="H2O")
        >>> participants = atp + water
        >>> len(participants.molecules)
        2
        >>> participants.molecules[0].name
        'ATP'
        >>> participants.molecules[1].name
        'H2O'
    """

    chebi_id: Optional[str] = None
    smiles: Optional[str] = None
    name: Optional[str] = None
    location: Optional[str] = None
    count: int = 1

    def __mul__(self, count: int) -> ParticipantBuilder:
        """Multiply participant by stoichiometric coefficient."""
        return ParticipantBuilder(
            chebi_id=self.chebi_id,
            smiles=self.smiles,
            name=self.name,
            location=self.location,
            count=self.count * count,
        )

    def __rmul__(self, count: int) -> ParticipantBuilder:
        """Right multiplication for stoichiometry (2 * H2O)."""
        return self.__mul__(count)

    def __getitem__(self, location: str) -> ParticipantBuilder:
        """Set location/compartment using subscript notation."""
        return ParticipantBuilder(
            chebi_id=self.chebi_id,
            smiles=self.smiles,
            name=self.name,
            location=location,
            count=self.count,
        )

    def __add__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ParticipantList:
        """Add participant builders together to form a participant list."""
        if isinstance(other, ParticipantBuilder):
            return ParticipantList([self, other])
        elif isinstance(other, ParticipantList):
            return ParticipantList([self] + other.molecules)
        else:
            raise TypeError(f"Cannot add ParticipantBuilder and {type(other)}")

    def __rshift__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a reaction using >> operator (left-to-right)."""
        left = ParticipantList([self])
        if isinstance(other, ParticipantBuilder):
            right = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            right = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(left, right, direction="left-to-right")

    def __lshift__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a reaction using << operator (right-to-left)."""
        right = ParticipantList([self])
        if isinstance(other, ParticipantBuilder):
            left = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            left = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(left, right, direction="right-to-left")

    def __or__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a bidirectional reaction using | operator."""
        left = ParticipantList([self])
        if isinstance(other, ParticipantBuilder):
            right = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            right = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(left, right, direction="bidirectional")

    def to_participant(self) -> Participant:
        """Convert to a Participant object."""
        return Participant(
            chebi_id=self.chebi_id,
            smiles=self.smiles,
            name=self.name,
            location=self.location,
            count=self.count,
        )


@dataclass
class ParticipantList:
    """A list of participant builders on one side of a reaction.

    This is created when participant builders are added together with the + operator.
    It's iterable and can be further extended with more participant builders.

    Examples:
        >>> from autarch.formula import ParticipantBuilder, ParticipantList
        >>>
        >>> # Create participant builders
        >>> atp = ParticipantBuilder(chebi_id="CHEBI:15422", name="ATP")
        >>> h2o = ParticipantBuilder(chebi_id="CHEBI:15377", name="H2O")
        >>> mg = ParticipantBuilder(chebi_id="CHEBI:18420", name="Mg2+")
        >>>
        >>> # Combine participant builders into a participant list
        >>> participants = atp + h2o
        >>> isinstance(participants, ParticipantList)
        True
        >>> len(participants)
        2
        >>>
        >>> # Add more participant builders
        >>> participants = participants + mg
        >>> len(participants)
        3
        >>>
        >>> # Iterate over participant builders
        >>> names = [mol.name for mol in participants]
        >>> names
        ['ATP', 'H2O', 'Mg2+']
        >>>
        >>> # Convert to Participant objects
        >>> participant_objs = participants.to_participants()
        >>> len(participant_objs)
        3
        >>> participant_objs[0].chebi_id
        'CHEBI:15422'
        >>>
        >>> # Use in reactions
        >>> adp = ParticipantBuilder(chebi_id="CHEBI:16761", name="ADP")
        >>> pi = ParticipantBuilder(chebi_id="CHEBI:43474", name="Pi")
        >>> rxn = participants >> (adp + pi + mg)
        >>> rxn.direction
        'left-to-right'
        >>> len(rxn.left)
        3
        >>> len(rxn.right)
        3
    """

    molecules: List[ParticipantBuilder] = field(default_factory=list)

    def __add__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ParticipantList:
        """Add more participant builders to this list."""
        if isinstance(other, ParticipantBuilder):
            return ParticipantList(self.molecules + [other])
        elif isinstance(other, ParticipantList):
            return ParticipantList(self.molecules + other.molecules)
        else:
            raise TypeError(f"Cannot add ParticipantList and {type(other)}")

    def __radd__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ParticipantList:
        """Right addition for participant builders."""
        if isinstance(other, ParticipantBuilder):
            return ParticipantList([other] + self.molecules)
        elif isinstance(other, ParticipantList):
            return ParticipantList(other.molecules + self.molecules)
        else:
            raise TypeError(f"Cannot add {type(other)} and ParticipantList")

    def __rshift__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a reaction using >> operator (left-to-right)."""
        if isinstance(other, ParticipantBuilder):
            right = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            right = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(self, right, direction="left-to-right")

    def __lshift__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a reaction using << operator (right-to-left)."""
        if isinstance(other, ParticipantBuilder):
            left = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            left = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(left, self, direction="right-to-left")

    def __or__(
        self, other: Union[ParticipantBuilder, ParticipantList]
    ) -> ReactionBuilder:
        """Create a bidirectional reaction using | operator."""
        if isinstance(other, ParticipantBuilder):
            right = ParticipantList([other])
        elif isinstance(other, ParticipantList):
            right = other
        else:
            raise TypeError(f"Cannot create reaction with {type(other)}")
        return ReactionBuilder(self, right, direction="bidirectional")

    def __iter__(self):
        """Iterate over participant builders in the list."""
        return iter(self.molecules)

    def __len__(self):
        """Number of distinct participant builders (not counting stoichiometry)."""
        return len(self.molecules)

    def to_participants(self) -> List[Participant]:
        """Convert all molecules to Participant objects."""
        return [mol.to_participant() for mol in self.molecules]


@dataclass
class ReactionBuilder:
    """Builder for creating Reaction objects with nice syntax.

    This is created by the >>, <<, or | operators and can be converted
    to a standard Reaction object for use with the classifier system.

    Examples:
        >>> from autarch.formula import ParticipantBuilder, ReactionBuilder
        >>>
        >>> # Create participant builders
        >>> atp = ParticipantBuilder(chebi_id="CHEBI:15422", name="ATP")
        >>> adp = ParticipantBuilder(chebi_id="CHEBI:16761", name="ADP")
        >>> h2o = ParticipantBuilder(chebi_id="CHEBI:15377", name="H2O")
        >>> pi = ParticipantBuilder(chebi_id="CHEBI:43474", name="Pi")
        >>>
        >>> # Create a reaction
        >>> rxn_builder = atp + h2o >> adp + pi
        >>> isinstance(rxn_builder, ReactionBuilder)
        True
        >>> rxn_builder.direction
        'left-to-right'
        >>>
        >>> # Convert to Reaction object
        >>> reaction = rxn_builder.build()
        >>> len(reaction.left_participants)
        2
        >>> len(reaction.right_participants)
        2
        >>>
        >>> # String representation
        >>> str(rxn_builder)
        'ATP + H2O → ADP + Pi'
        >>>
        >>> # Bidirectional reaction
        >>> nad = ParticipantBuilder(chebi_id="CHEBI:57540", name="NAD+")
        >>> nadh = ParticipantBuilder(chebi_id="CHEBI:57945", name="NADH")
        >>> h_plus = ParticipantBuilder(chebi_id="CHEBI:15378", name="H+")
        >>> redox = nad + h_plus | nadh
        >>> redox.direction
        'bidirectional'
        >>> str(redox)
        'NAD+ + H+ ⇌ NADH'
        >>>
        >>> # Reverse reaction
        >>> reverse = adp + pi << atp + h2o
        >>> reverse.direction
        'right-to-left'
        >>> str(reverse)
        'ATP + H2O ← ADP + Pi'
    """

    left: ParticipantList
    right: ParticipantList
    direction: str = "bidirectional"

    def build(self) -> Reaction:
        """Build a Reaction object from this builder."""
        return Reaction(
            left_participants=self.left.to_participants(),
            right_participants=self.right.to_participants(),
        )

    def __str__(self) -> str:
        """String representation of the reaction."""
        arrow = {"left-to-right": "→", "right-to-left": "←", "bidirectional": "⇌"}.get(
            self.direction, "→"
        )

        left_str = " + ".join(
            f"{m.count} {m.name or m.chebi_id or 'molecule'}"
            if m.count > 1
            else f"{m.name or m.chebi_id or 'molecule'}"
            for m in self.left.molecules
        )
        right_str = " + ".join(
            f"{m.count} {m.name or m.chebi_id or 'molecule'}"
            if m.count > 1
            else f"{m.name or m.chebi_id or 'molecule'}"
            for m in self.right.molecules
        )

        return f"{left_str} {arrow} {right_str}"


# Factory function for creating participant builders
def molecule(
    chebi_id: Optional[str] = None,
    smiles: Optional[str] = None,
    name: Optional[str] = None,
) -> ParticipantBuilder:
    """Create a participant builder for use in reaction formulas.

    This is a convenience factory function for creating ParticipantBuilder objects
    with a simpler syntax than using the class constructor directly.

    Args:
        chebi_id: CHEBI identifier (e.g., "CHEBI:15422")
        smiles: SMILES string representation
        name: Common name for the molecule

    Returns:
        ParticipantBuilder object that can be used with operators

    Examples:
        >>> # Create with CHEBI ID and name
        >>> atp = molecule("CHEBI:15422", name="ATP")
        >>> atp.chebi_id
        'CHEBI:15422'
        >>> atp.name
        'ATP'
        >>>
        >>> # Create with SMILES
        >>> ethanol = molecule(smiles="CCO", name="ethanol")
        >>> ethanol.smiles
        'CCO'
        >>> ethanol.name
        'ethanol'
        >>>
        >>> # Use in a reaction
        >>> adp = molecule("CHEBI:16761", name="ADP")
        >>> water = molecule("CHEBI:15377", name="H2O")
        >>> pi = molecule("CHEBI:43474", name="Pi")
        >>> hydrolysis = atp + water >> adp + pi
        >>> str(hydrolysis)
        'ATP + H2O → ADP + Pi'
        >>>
        >>> # Can be used with all operators
        >>> atp2 = atp * 2
        >>> atp2.count
        2
        >>> atp_cytoplasm = atp["cytoplasm"]
        >>> atp_cytoplasm.location
        'cytoplasm'
    """
    return ParticipantBuilder(chebi_id=chebi_id, smiles=smiles, name=name)


# Convenience function for creating reactions
def rxn(formula_string: str) -> Reaction:
    """Parse a reaction formula string (future enhancement).

    This could parse strings like "ATP + H2O -> ADP + Pi"
    but for now just raises NotImplementedError.
    """
    raise NotImplementedError(
        "String parsing not yet implemented. Use operator syntax instead."
    )


# Pre-defined common participant builders for backward compatibility
# Note: These are now also available in autarch.vocabulary with many more molecules
ATP = molecule("CHEBI:15422", name="ATP")
ADP = molecule("CHEBI:16761", name="ADP")
AMP = molecule("CHEBI:16027", name="AMP")
Pi = molecule("CHEBI:43474", name="Pi")  # Inorganic phosphate
PPi = molecule("CHEBI:33019", name="PPi")  # Pyrophosphate
H2O = molecule("CHEBI:15377", name="H2O")
H = molecule("CHEBI:15378", name="H+")
NAD = molecule("CHEBI:57540", name="NAD+")
NADH = molecule("CHEBI:57945", name="NADH")
NADP = molecule("CHEBI:58349", name="NADP+")
NADPH = molecule("CHEBI:57783", name="NADPH")
FAD = molecule("CHEBI:57692", name="FAD")
FADH2 = molecule("CHEBI:57618", name="FADH2")
CoA = molecule("CHEBI:15346", name="CoA")
glucose = molecule("CHEBI:17234", name="glucose")
O2 = molecule("CHEBI:15379", name="O2")
CO2 = molecule("CHEBI:16526", name="CO2")


def demonstrate_classification():
    """Demonstrate how the formula DSL integrates with the classifier system.

    This function shows how reactions created with the formula DSL can be
    classified using the ReactionClassifier.

    Examples:
        >>> from autarch.formula import ATP, ADP, Pi, H2O, NAD, NADH, H
        >>> from autarch.classifier import ReactionClassifier
        >>>
        >>> # Create a hydrolase reaction (ATP hydrolysis)
        >>> hydrolysis = ATP + H2O >> ADP + Pi
        >>> reaction = hydrolysis.build()
        >>>
        >>> # Classify the reaction
        >>> classifier = ReactionClassifier()
        >>> results = classifier.classify(reaction)
        >>>
        >>> # Check if it's a hydrolase
        >>> 'Hydrolase' in results
        True
        >>> results['Hydrolase'].is_member
        True
        >>>
        >>> # Create an oxidoreductase reaction
        >>> oxidation = NAD + H >> NADH
        >>> reaction = oxidation.build()
        >>> results = classifier.classify(reaction)
        >>>
        >>> # Check classifications
        >>> 'Oxidoreductase' in results
        True
        >>> results['Oxidoreductase'].is_member
        True
        >>>
        >>> # Transport reaction with locations
        >>> transport = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
        >>> reaction = transport.build()
        >>>
        >>> # Check it's recognized as transport
        >>> reaction.is_transport_reaction()
        True
        >>>
        >>> # Get transport details
        >>> transported = reaction.get_transported_molecules()
        >>> len(transported)
        2
        >>>
        >>> # Complex reaction with stoichiometry
        >>> complex_rxn = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)
        >>> reaction = complex_rxn.build()
        >>>
        >>> # Check participant counts
        >>> reaction.left_participants[1].count  # O2
        6
        >>> reaction.right_participants[0].count  # CO2
        6
        >>>
        >>> # Print a reaction nicely
        >>> print(hydrolysis)
        ATP + H2O → ADP + Pi
        >>> print(oxidation)
        NAD+ + H+ → NADH
        >>> print(transport)
        ATP + ADP → ATP + ADP
    """
    pass  # Function exists just for its doctest
