"""Data models for autarch.

This module defines the core data structures used throughout the autarch
reaction classification system. It includes models for chemical participants,
reactions, and classification results.

The data model supports:
- Chemical species with CHEBI IDs, SMILES, and stoichiometry
- Polymer types for reactions that can't be represented as SMILES
- Cellular compartments/locations for transport reactions
- Reaction representation with left and right participants
- RHEA and GO term metadata
- Classification results with explanations
"""

from enum import Enum
from typing import TYPE_CHECKING, Optional
from pydantic import BaseModel, Field
from rdkit import Chem

if TYPE_CHECKING:
    from autarch.moiety import Moiety


class PolymerType(str, Enum):
    """Types of biological polymers in reactions.

    Used for reactions involving polymeric substrates that cannot be
    fully represented as SMILES strings. The polymer_index field on
    Participant indicates the length (n, n+1, n-1).

    Examples:
        >>> PolymerType.RNA.value
        'rna'
        >>> PolymerType.GLUCAN in [PolymerType.GLUCAN, PolymerType.CHITIN]
        True
    """

    # Nucleic acids
    RNA = "rna"
    DNA = "dna"
    TRNA = "trna"
    MRNA = "mrna"
    RRNA = "rrna"

    # Proteins
    PROTEIN = "protein"
    PEPTIDE = "peptide"

    # Polysaccharides - glucose-based
    GLUCAN = "glucan"  # Generic glucose polymer (starch, glycogen, cellulose)
    STARCH = "starch"
    GLYCOGEN = "glycogen"
    CELLULOSE = "cellulose"

    # Polysaccharides - other
    CHITIN = "chitin"  # N-acetyl-D-glucosamine polymer
    CHITOSAN = "chitosan"  # Deacetylated chitin
    PEPTIDOGLYCAN = "peptidoglycan"  # Bacterial cell wall
    SIALIC_ACID_POLYMER = "sialic_acid_polymer"  # Polysialic acid
    GALACTURONAN = "galacturonan"  # Pectin backbone
    FRUCTAN = "fructan"  # Fructose polymer (inulin, levan)

    # Other polymers
    PHOSPHATE_CHAIN = "phosphate_chain"  # Polyphosphate
    TEICHOIC_ACID = "teichoic_acid"  # Bacterial cell wall component
    POLYSULFUR = "polysulfur"  # Sulfur chain

    # Generic fallbacks
    POLYSACCHARIDE = "polysaccharide"
    POLYNUCLEOTIDE = "polynucleotide"
    POLYPEPTIDE = "polypeptide"
    OTHER = "other"


class ClassificationResult(BaseModel):
    """Result of checking if a reaction belongs to a class.

    This simple class encapsulates the binary result of a classification
    check along with a human-readable explanation of the reasoning.

    Attributes:
        is_member: Whether the reaction belongs to this class
        explanation: Human-readable explanation of the classification decision

    Examples:
        >>> result = ClassificationResult(
        ...     is_member=True,
        ...     explanation="Hydrolysis: water consumed, fragmentation detected"
        ... )
        >>> result.is_member
        True
        >>> result.explanation
        'Hydrolysis: water consumed, fragmentation detected'
    """

    is_member: bool
    explanation: str


class Participant(BaseModel):
    """A chemical species participating in a reaction.

    Represents a molecule/ion with its structure, identifier, stoichiometry,
    and optionally its cellular location/compartment.

    Attributes:
        chebi_id: CHEBI identifier (e.g., "CHEBI:15377" for water)
        smiles: SMILES string representation of the molecule
        count: Stoichiometric coefficient (default 1)
        stoichiometry: Variable stoichiometry expression (e.g., "n", "n+1", "2n")
            Used for polymer reactions where stoichiometry depends on polymer length.
            When set, this takes precedence over count for display purposes.
        polymer_index: Polymer length indicator (e.g., "n", "n+1", "n-1")
            Used when the participant itself is a polymer whose length changes.
            Example: tRNA(n) -> tRNA(n+1) or [(glucosyl)](n) -> [(glucosyl)](n+1)
        location: Cellular compartment or location (e.g., "cytoplasm", "mitochondria", "in", "out")
        name: Common name of the molecule (optional)
        formula: Molecular formula (optional, e.g., "H2O")
        charge: Net charge of the species (optional)
        mol: RDKit Mol object (excluded from serialization)

    Examples:
        >>> # Regular participant with fixed stoichiometry
        >>> atp = Participant(name="ATP", count=2)
        >>> str(atp)
        '2 ATP'
        >>>
        >>> # Participant with variable stoichiometry (polymer reaction)
        >>> atp_n = Participant(name="ATP", stoichiometry="n")
        >>> str(atp_n)
        'n ATP'
        >>>
        >>> # Polymer that grows
        >>> trna = Participant(name="tRNA", polymer_index="n+1")
        >>> str(trna)
        'tRNA(n+1)'
    """

    # Core identifiers
    chebi_id: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None  # InChI for structural comparison (stereochemistry)

    # Stoichiometry and properties
    count: int = 1
    stoichiometry: Optional[str] = None  # Variable stoichiometry: "n", "2n", "n+1"
    polymer_index: Optional[str] = None  # Polymer length: "n", "n+1", "n-1"
    location: Optional[str] = None
    name: Optional[str] = None
    formula: Optional[str] = None
    charge: Optional[int] = None

    # Polymer information (for reactions that can't be represented as SMILES)
    polymer_type: Optional[PolymerType] = None  # Type of polymer (RNA, GLUCAN, etc.)
    linkage: Optional[str] = None  # Glycosidic/phosphodiester linkage (e.g., "1->4", "2->8")
    monomer: Optional[str] = None  # Repeating unit description

    # RDKit Mol object (not serializable)
    mol: Optional[Chem.Mol] = Field(default=None, exclude=True)

    model_config = {"arbitrary_types_allowed": True}  # Allow RDKit Mol objects

    def get_mol(self) -> Optional[Chem.Mol]:
        """Get or create RDKit Mol object.

        Lazily creates an RDKit Mol object from the SMILES string if not already
        created. This allows for chemical structure analysis while avoiding the
        overhead of creating Mol objects when not needed.

        Returns:
            RDKit Mol object if SMILES is available, None otherwise

        Examples:
            >>> p = Participant(smiles="CCO", name="ethanol")
            >>> mol = p.get_mol()
            >>> mol.GetNumAtoms()
            3
        """
        if self.mol is None and self.smiles:
            self.mol = Chem.MolFromSmiles(self.smiles)
        return self.mol

    def is_same_molecule(self, other: "Participant") -> bool:
        """Check if this participant represents the same molecule as another.

        Two participants are considered the same molecule if they have the same
        CHEBI ID or the same SMILES string. Location/compartment information
        is ignored for this comparison, as the same molecule in different
        locations is still chemically identical.

        Args:
            other: Another Participant to compare with

        Returns:
            True if the participants represent the same molecule

        Examples:
            >>> atp1 = Participant(chebi_id="CHEBI:15422", location="cytoplasm")
            >>> atp2 = Participant(chebi_id="CHEBI:15422", location="mitochondria")
            >>> atp1.is_same_molecule(atp2)
            True
            >>>
            >>> water = Participant(chebi_id="CHEBI:15377")
            >>> atp1.is_same_molecule(water)
            False
        """
        if self.chebi_id and other.chebi_id:
            return self.chebi_id == other.chebi_id
        if self.smiles and other.smiles:
            return self.smiles == other.smiles
        return False

    def has_location(self) -> bool:
        """Check if this participant has a specified location.

        Returns:
            True if location is not None

        Examples:
            >>> p1 = Participant(chebi_id="CHEBI:15422", location="cytoplasm")
            >>> p1.has_location()
            True
            >>> p2 = Participant(chebi_id="CHEBI:15422")
            >>> p2.has_location()
            False
        """
        return self.location is not None

    def has_variable_stoichiometry(self) -> bool:
        """Check if this participant has variable (n-based) stoichiometry.

        Returns:
            True if stoichiometry contains 'n' (e.g., "n", "2n", "n+1")

        Examples:
            >>> p1 = Participant(name="ATP", stoichiometry="n")
            >>> p1.has_variable_stoichiometry()
            True
            >>> p2 = Participant(name="ATP", count=2)
            >>> p2.has_variable_stoichiometry()
            False
        """
        return self.stoichiometry is not None and "n" in self.stoichiometry

    def is_polymer(self) -> bool:
        """Check if this participant is a polymer with variable length.

        Returns:
            True if polymer_index is set (e.g., "n", "n+1", "n-1")

        Examples:
            >>> p1 = Participant(name="tRNA", polymer_index="n+1")
            >>> p1.is_polymer()
            True
            >>> p2 = Participant(name="ATP")
            >>> p2.is_polymer()
            False
        """
        return self.polymer_index is not None

    # Moiety detection methods (SMARTS-based)
    def has_moiety(self, moiety: "Moiety") -> bool:
        """Check if this participant contains a moiety using SMARTS.

        Args:
            moiety: Moiety enum value to check for

        Returns:
            True if the moiety is present in this molecule

        Examples:
            >>> from autarch.moiety import Moiety
            >>> p = Participant(smiles="CC(=O)SC", chebi_id="test")  # thioester
            >>> p.has_moiety(Moiety.THIOESTER)
            True
        """
        if not self.smiles:
            return False
        from autarch.moiety import has_moiety as _has_moiety
        return _has_moiety(self.smiles, moiety, self.chebi_id)

    def is_phosphorylated(self) -> bool:
        """Check if this participant contains a phosphate group.

        Examples:
            >>> p = Participant(smiles="OP(O)(O)=O")  # phosphoric acid
            >>> p.is_phosphorylated()
            True
            >>> p2 = Participant(smiles="O")  # water
            >>> p2.is_phosphorylated()
            False
        """
        if not self.smiles:
            return False
        from autarch.moiety import is_phosphorylated as _is_phosphorylated
        return _is_phosphorylated(self.smiles, self.chebi_id)

    def is_thioester(self) -> bool:
        """Check if this participant contains a thioester bond (CoA-like).

        Detects C(=O)-S linkages like those found in acetyl-CoA,
        palmitoyl-CoA, and other acyl-CoA derivatives.

        Examples:
            >>> p = Participant(smiles="CC(=O)SC")  # methyl thioacetate
            >>> p.is_thioester()
            True
            >>> p2 = Participant(smiles="CC(=O)OC")  # methyl acetate (ester)
            >>> p2.is_thioester()
            False
        """
        if not self.smiles:
            return False
        from autarch.moiety import is_thioester as _is_thioester
        return _is_thioester(self.smiles, self.chebi_id)

    def is_glycoside(self) -> bool:
        """Check if this participant contains a glycosidic bond.

        Detects both O-glycosides and N-glycosides.
        """
        if not self.smiles:
            return False
        from autarch.moiety import is_glycoside as _is_glycoside
        return _is_glycoside(self.smiles, self.chebi_id)

    def is_sulfated(self) -> bool:
        """Check if this participant contains a sulfate group."""
        if not self.smiles:
            return False
        from autarch.moiety import is_sulfated as _is_sulfated
        return _is_sulfated(self.smiles, self.chebi_id)

    def is_sulfamate(self) -> bool:
        """Check if this participant contains a sulfamate group (S-N bond).

        Sulfamates have a sulfur-nitrogen bond, characteristic of compounds
        like cyclohexylsulfamate or N-sulfo-D-glucosamine.
        """
        if not self.smiles:
            return False
        from autarch.moiety import is_sulfamate as _is_sulfamate
        return _is_sulfamate(self.smiles, self.chebi_id)

    def is_carboxylic_acid(self) -> bool:
        """Check if this participant contains a carboxylic acid group."""
        if not self.smiles:
            return False
        from autarch.moiety import is_carboxylic_acid as _is_carboxylic_acid
        return _is_carboxylic_acid(self.smiles, self.chebi_id)

    def __str__(self) -> str:
        """String representation of the participant.

        Creates a human-readable string showing the stoichiometry, name/ID,
        polymer index, and location of the participant.

        Returns:
            String representation

        Examples:
            >>> str(Participant(name="ATP", count=2, location="cytoplasm"))
            '2 ATP (cytoplasm)'
            >>> str(Participant(chebi_id="CHEBI:15422"))
            'CHEBI:15422'
            >>> str(Participant(name="ATP", stoichiometry="n"))
            'n ATP'
            >>> str(Participant(name="tRNA", polymer_index="n+1"))
            'tRNA(n+1)'
        """
        parts = []

        # Stoichiometry: variable (n, 2n) takes precedence over count
        if self.stoichiometry:
            parts.append(self.stoichiometry)
        elif self.count > 1:
            parts.append(f"{self.count}")

        # Name/identifier with optional polymer index
        name_part = None
        if self.name:
            name_part = self.name
        elif self.chebi_id:
            name_part = self.chebi_id
        elif self.smiles:
            name_part = self.smiles[:20] + "..." if len(self.smiles) > 20 else self.smiles

        if name_part:
            if self.polymer_index:
                name_part = f"{name_part}({self.polymer_index})"
            parts.append(name_part)

        # Location
        if self.location:
            parts.append(f"({self.location})")

        return " ".join(parts) if parts else "Unknown participant"


class Reaction(BaseModel):
    """A chemical reaction with left and right participants.

    Represents a balanced chemical equation with reactants (left) and products (right).
    Can include location information for transport reactions. This is the core
    data structure used throughout the classification system.

    Attributes:
        left_participants: List of reactants (left side of equation)
        right_participants: List of products (right side of equation)

    Examples:
        >>> from autarch.datamodel import Reaction, Participant
        >>>
        >>> # ATP hydrolysis
        >>> reaction = Reaction(
        ...     left_participants=[
        ...         Participant(chebi_id="CHEBI:15422", name="ATP"),
        ...         Participant(chebi_id="CHEBI:15377", name="water")
        ...     ],
        ...     right_participants=[
        ...         Participant(chebi_id="CHEBI:16761", name="ADP"),
        ...         Participant(chebi_id="CHEBI:43474", name="phosphate")
        ...     ]
        ... )
        >>> len(reaction.left_participants)
        2
        >>> len(reaction.right_participants)
        2
        >>>
        >>> # Transport reaction
        >>> transport = Reaction(
        ...     left_participants=[
        ...         Participant(chebi_id="CHEBI:15422", location="out"),
        ...         Participant(chebi_id="CHEBI:16761", location="in")
        ...     ],
        ...     right_participants=[
        ...         Participant(chebi_id="CHEBI:15422", location="in"),
        ...         Participant(chebi_id="CHEBI:16761", location="out")
        ...     ]
        ... )
        >>> transport.is_transport_reaction()
        True
    """

    left_participants: list[Participant] = Field(default_factory=list)
    right_participants: list[Participant] = Field(default_factory=list)
    label: str = ""  # Human-readable equation string from RHEA

    def all_participants(self) -> list[Participant]:
        """Get all participants from both sides.

        Returns:
            Combined list of all participants from left and right sides

        Examples:
            >>> reaction = Reaction(
            ...     left_participants=[Participant(name="A"), Participant(name="B")],
            ...     right_participants=[Participant(name="C")]
            ... )
            >>> len(reaction.all_participants())
            3
        """
        return self.left_participants + self.right_participants

    def is_transport_reaction(self) -> bool:
        """Check if this is a transport reaction.

        A transport reaction has the same molecules on both sides but in different
        locations/compartments. This is detected by comparing molecules while
        ignoring their locations and checking if locations differ.

        Returns:
            True if this is a transport reaction

        Examples:
            >>> # ATP/ADP antiporter
            >>> transport = Reaction(
            ...     left_participants=[
            ...         Participant(chebi_id="CHEBI:15422", location="out"),
            ...         Participant(chebi_id="CHEBI:16761", location="in")
            ...     ],
            ...     right_participants=[
            ...         Participant(chebi_id="CHEBI:15422", location="in"),
            ...         Participant(chebi_id="CHEBI:16761", location="out")
            ...     ]
            ... )
            >>> transport.is_transport_reaction()
            True
            >>>
            >>> # Regular reaction (not transport)
            >>> regular = Reaction(
            ...     left_participants=[
            ...         Participant(chebi_id="CHEBI:15422"),
            ...         Participant(chebi_id="CHEBI:15377")
            ...     ],
            ...     right_participants=[
            ...         Participant(chebi_id="CHEBI:16761"),
            ...         Participant(chebi_id="CHEBI:43474")
            ...     ]
            ... )
            >>> regular.is_transport_reaction()
            False
        """
        # Check if any participants have locations
        if not any(p.has_location() for p in self.all_participants()):
            return False

        # Group participants by molecule (ignoring location)
        left_molecules: dict[tuple[Optional[str], int], list[Optional[str]]] = {}
        right_molecules: dict[tuple[Optional[str], int], list[Optional[str]]] = {}

        for p in self.left_participants:
            key = (p.chebi_id or p.smiles, p.count)
            if key not in left_molecules:
                left_molecules[key] = []
            left_molecules[key].append(p.location)

        for p in self.right_participants:
            key = (p.chebi_id or p.smiles, p.count)
            if key not in right_molecules:
                right_molecules[key] = []
            right_molecules[key].append(p.location)

        # Check if same molecules appear on both sides with different locations
        for key in left_molecules:
            if key in right_molecules:
                left_locs = set(left_molecules[key])
                right_locs = set(right_molecules[key])
                if left_locs != right_locs:
                    return True

        return False

    def get_transported_molecules(self) -> list[tuple[str, str, str]]:
        """Get list of transported molecules with their locations.

        Identifies molecules that appear on both sides of the reaction with
        different locations, indicating transport across compartments.

        Returns:
            List of tuples: (molecule_id, from_location, to_location)

        Examples:
            >>> transport = Reaction(
            ...     left_participants=[
            ...         Participant(chebi_id="CHEBI:15422", name="ATP", location="out"),
            ...         Participant(chebi_id="CHEBI:16761", name="ADP", location="in")
            ...     ],
            ...     right_participants=[
            ...         Participant(chebi_id="CHEBI:15422", name="ATP", location="in"),
            ...         Participant(chebi_id="CHEBI:16761", name="ADP", location="out")
            ...     ]
            ... )
            >>> transported_mols = transport.get_transported_molecules()
            >>> len(transported_mols)
            2
            >>> # ATP moves from out to in
            >>> transported_mols[0]
            ('CHEBI:15422', 'out', 'in')
            >>> # ADP moves from in to out
            >>> transported_mols[1]
            ('CHEBI:16761', 'in', 'out')
        """
        transported = []

        for left_p in self.left_participants:
            for right_p in self.right_participants:
                if (
                    left_p.is_same_molecule(right_p)
                    and left_p.location != right_p.location
                ):
                    mol_id = (
                        left_p.chebi_id or left_p.smiles or left_p.name or "unknown"
                    )
                    transported.append(
                        (
                            mol_id,
                            left_p.location or "unknown",
                            right_p.location or "unknown",
                        )
                    )

        return transported

    def has_participant_with_chebi(self, chebi_id: str) -> bool:
        """Check if reaction has a participant with given CHEBI ID.

        Args:
            chebi_id: CHEBI identifier to search for

        Returns:
            True if any participant has the given CHEBI ID

        Examples:
            >>> reaction = Reaction(
            ...     left_participants=[Participant(chebi_id="CHEBI:15422")],
            ...     right_participants=[Participant(chebi_id="CHEBI:16761")]
            ... )
            >>> reaction.has_participant_with_chebi("CHEBI:15422")
            True
            >>> reaction.has_participant_with_chebi("CHEBI:99999")
            False
        """
        return any(p.chebi_id == chebi_id for p in self.all_participants())

    def has_participant_with_name(self, name: str) -> bool:
        """Check if reaction has a participant with given name.

        Args:
            name: Name to search for

        Returns:
            True if any participant has the given name

        Examples:
            >>> reaction = Reaction(
            ...     left_participants=[Participant(name="ATP"), Participant(name="water")],
            ...     right_participants=[Participant(name="ADP"), Participant(name="phosphate")]
            ... )
            >>> reaction.has_participant_with_name("ATP")
            True
            >>> reaction.has_participant_with_name("glucose")
            False
        """
        return any(p.name == name for p in self.all_participants())

    def has_participant_with_smiles(self, smiles: str) -> bool:
        """Check if reaction has a participant with given SMILES.

        Args:
            smiles: SMILES string to search for

        Returns:
            True if any participant has the given SMILES

        Examples:
            >>> reaction = Reaction(
            ...     left_participants=[Participant(smiles="O"), Participant(smiles="CCO")],
            ...     right_participants=[Participant(smiles="CC=O")]
            ... )
            >>> reaction.has_participant_with_smiles("O")
            True
            >>> reaction.has_participant_with_smiles("CC")
            False
        """
        return any(p.smiles == smiles for p in self.all_participants())

    def __str__(self) -> str:
        """String representation of the reaction.

        Creates a human-readable equation showing all participants.

        Returns:
            String representation with arrow

        Examples:
            >>> reaction = Reaction(
            ...     left_participants=[
            ...         Participant(name="ATP"),
            ...         Participant(name="H2O")
            ...     ],
            ...     right_participants=[
            ...         Participant(name="ADP"),
            ...         Participant(name="Pi")
            ...     ]
            ... )
            >>> str(reaction)
            'ATP + H2O → ADP + Pi'
        """
        left = " + ".join(str(p) for p in self.left_participants)
        right = " + ".join(str(p) for p in self.right_participants)
        return f"{left} → {right}"


class RheaTerm(BaseModel):
    """Represents a RHEA reaction term with its metadata and reaction data.

    RHEA is a comprehensive resource of expert-curated biochemical reactions.
    This class encapsulates a RHEA reaction with all its associated metadata
    including EC numbers, GO terms, and directionality information.

    Attributes:
        rhea_id: RHEA reaction identifier (e.g., "RHEA:18777")
        label: Human-readable label for the reaction
        reaction: The actual Reaction object with participants
        ec_numbers: List of EC numbers that catalyze this reaction
        go_terms: List of GO molecular function terms
        direction: Reaction direction (bidirectional, left-to-right, right-to-left)
        parent_rhea_id: Parent reaction ID for directional variants

    Examples:
        >>> rhea_term = RheaTerm(
        ...     rhea_id="RHEA:15421",
        ...     label="ATP + H2O = ADP + phosphate",
        ...     direction="left-to-right",
        ...     ec_numbers=["3.6.1.3"],
        ...     go_terms=["GO:0016787"]
        ... )
        >>> rhea_term.rhea_id
        'RHEA:15421'
        >>> rhea_term.ec_numbers
        ['3.6.1.3']
        >>> rhea_term.direction
        'left-to-right'
    """

    rhea_id: str
    label: str = ""
    reaction: Optional[Reaction] = None
    ec_numbers: list[str] = Field(default_factory=list)
    go_terms: list[str] = Field(default_factory=list)
    direction: str = "bidirectional"
    parent_rhea_id: Optional[str] = None


class GoTerm(BaseModel):
    """Represents a GO term with its mappings to various databases.

    Gene Ontology (GO) terms describe molecular functions, biological processes,
    and cellular components. This class focuses on molecular function terms
    that are associated with enzymatic reactions.

    Attributes:
        go_id: GO term identifier (e.g., "GO:0004553")
        label: Human-readable label for the GO term
        definition: GO term definition
        rhea_ids: List of RHEA reaction identifiers
        ec_numbers: List of EC numbers
        ancestors: List of ancestor GO term IDs (closure)

    Examples:
        >>> go_term = GoTerm(
        ...     go_id="GO:0016787",
        ...     label="hydrolase activity",
        ...     definition="Catalysis of the hydrolysis of various bonds",
        ...     ec_numbers=["3.-.-.-"],
        ...     rhea_ids=["RHEA:15421", "RHEA:15422", "RHEA:15423"],
        ...     ancestors=["GO:0003824", "GO:0008150"]
        ... )
        >>> go_term.go_id
        'GO:0016787'
        >>> go_term.label
        'hydrolase activity'
        >>> len(go_term.rhea_ids)
        3
        >>> "3.-.-.-" in go_term.ec_numbers
        True
    """

    go_id: str
    label: str = ""
    definition: str = ""
    rhea_ids: list[str] = Field(default_factory=list)
    ec_numbers: list[str] = Field(default_factory=list)
    ancestors: list[str] = Field(default_factory=list)
