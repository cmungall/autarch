# Chemical Reaction Formula DSL

The autarch library provides a powerful Domain Specific Language (DSL) for expressing chemical reactions using Python's operator overloading. This provides a natural, readable syntax that closely mirrors how chemists write reactions.

## Quick Start

```python
from autarch.formula import ATP, ADP, Pi, H2O, H
from autarch.classifier import ReactionClassifier

# Express a reaction naturally
hydrolysis = ATP + H2O >> ADP + Pi

# Convert to Reaction object and classify
reaction = hydrolysis.build()
classifier = ReactionClassifier()
results = classifier.classify(reaction)

print(f"Is hydrolase? {results['Hydrolase'].is_member}")
# Output: Is hydrolase? True
```

## Core Operators

### Combining Molecules: `+`
```python
# Combine molecules to form reaction participants
reactants = ATP + H2O + Mg
products = ADP + Pi + Mg
```

### Stoichiometry: `*`
```python
# Specify stoichiometric coefficients
glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)

# Both forms work
2 * H2O  # or  H2O * 2
```

### Compartments/Locations: `[]`
```python
# Specify cellular compartments
ATP["cytoplasm"] + substrate >> ADP["cytoplasm"] + product

# Transport reactions
(ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
```

### Reaction Directions

- `>>` : Left-to-right (forward)
- `<<` : Right-to-left (reverse)  
- `|` : Bidirectional (equilibrium)

```python
# Forward reaction
ATP + H2O >> ADP + Pi

# Reverse reaction (note: this reverses the actual direction)
ADP + Pi << ATP + H2O

# Bidirectional/equilibrium
NAD + H | NADH
```

## Complete Examples

### 1. ATP Hydrolysis (Hydrolase)
```python
from autarch.formula import ATP, ADP, Pi, H2O

reaction = (ATP + H2O >> ADP + Pi).build()
# This is recognized as a hydrolase reaction
```

### 2. Cellular Respiration (with stoichiometry)
```python
from autarch.formula import glucose, O2, CO2, H2O

respiration = glucose + (O2 * 6) >> (CO2 * 6) + (H2O * 6)
reaction = respiration.build()
```

### 3. NAD/NADH Redox (Oxidoreductase)
```python
from autarch.formula import NAD, NADH, H

# Bidirectional redox reaction
redox = NAD + H | NADH
reaction = redox.build()
```

### 4. ATP/ADP Translocase (Transport)
```python
from autarch.formula import ATP, ADP

# Antiporter - ATP goes in, ADP goes out
transport = (ATP["out"] + ADP["in"]) >> (ATP["in"] + ADP["out"])
reaction = transport.build()

# Check if it's a transport reaction
print(reaction.is_transport_reaction())  # True

# Get transport details
for mol_id, from_loc, to_loc in reaction.get_transported_molecules():
    print(f"{mol_id}: {from_loc} → {to_loc}")
```

### 5. Kinase Reaction (Phosphorylation)
```python
from autarch.formula import ATP, ADP, molecule

# Create custom substrate
substrate = molecule("CHEBI:17234", name="glucose")
product = molecule("CHEBI:17665", name="glucose-6-P")

# Phosphorylation reaction
kinase = ATP + substrate >> ADP + product
reaction = kinase.build()
```

### 6. Custom Molecules
```python
from autarch.formula import molecule, NAD, NADH, H

# Create molecules from CHEBI IDs
ethanol = molecule("CHEBI:16236", name="ethanol")
acetaldehyde = molecule("CHEBI:15343", name="acetaldehyde")

# Or from SMILES
ethanol = molecule(smiles="CCO", name="ethanol")
acetaldehyde = molecule(smiles="CC=O", name="acetaldehyde")

# Use in reactions
alcohol_dehydrogenase = ethanol + NAD >> acetaldehyde + NADH + H
```

## Integration with Classifiers

The DSL seamlessly integrates with the reaction classification system:

```python
from autarch.formula import ATP, ADP, Pi, H2O
from autarch.classifier import ReactionClassifier

# Create reaction using DSL
reaction_formula = ATP + H2O >> ADP + Pi

# Build to Reaction object
reaction = reaction_formula.build()

# Classify
classifier = ReactionClassifier()
results = classifier.classify(reaction)

# Check specific classes
if results['Hydrolase'].is_member:
    print("This is a hydrolase reaction")
    print(f"Reason: {results['Hydrolase'].explanation}")

# Get all matching classes
matching = classifier.get_matching_classes(reaction)
print(f"Matches: {', '.join(matching)}")
```

## Pre-defined Molecules

The library includes common molecules:

- **Energy**: ATP, ADP, AMP, Pi, PPi
- **Water/Protons**: H2O, H (H+)
- **NAD System**: NAD (NAD+), NADH, NADP (NADP+), NADPH
- **FAD System**: FAD, FADH2
- **Others**: CoA, glucose, O2, CO2

## Advanced Usage

### Complex Metabolic Pathways
```python
# Glycolysis steps
step1 = glucose + ATP >> glucose_6P + ADP + H
step2 = glucose_6P >> fructose_6P
step3 = fructose_6P + ATP >> fructose_16BP + ADP + H

# TCA cycle reactions
citrate_synthase = acetyl_CoA + oxaloacetate + H2O >> citrate + CoA + H
isocitrate_dh = isocitrate + NAD >> alpha_KG + NADH + CO2
```

### Operator Precedence
```python
# Multiplication binds tighter than addition
ATP * 2 + H2O  # Correct: 2 ATP molecules plus 1 H2O
# NOT: (ATP + H2O) * 2  # This would fail - can't multiply ParticipantList
```

### String Representation
```python
rxn = (ATP * 2) + H2O >> (ADP * 2) + (Pi * 2)
print(rxn)
# Output: 2 ATP + H2O → 2 ADP + 2 Pi

rxn = NAD + H | NADH
print(rxn)
# Output: NAD+ + H+ ⇌ NADH
```

## Benefits of the DSL

1. **Readable**: Reactions look like chemical equations
2. **Type-safe**: Full IDE support with autocomplete
3. **Validated**: Integrated with the classification system
4. **Flexible**: Supports stoichiometry, compartments, and directions
5. **Extensible**: Easy to add custom molecules

## Migration from Manual Construction

Instead of:
```python
# Old way - manual construction
reaction = Reaction(
    left_participants=[
        Participant(chebi_id="CHEBI:15422", name="ATP"),
        Participant(chebi_id="CHEBI:15377", name="H2O")
    ],
    right_participants=[
        Participant(chebi_id="CHEBI:16761", name="ADP"),
        Participant(chebi_id="CHEBI:43474", name="Pi")
    ]
)
```

Use:
```python
# New way - DSL
reaction = (ATP + H2O >> ADP + Pi).build()
```

The DSL provides a more natural and maintainable way to express chemical reactions while maintaining full compatibility with the existing autarch infrastructure.