# Fundamental Challenges in Reaction Classification

This document reflects on the goals, assumptions, and limitations of using declarative structural rules to classify enzyme-catalyzed reactions.

## The Original Vision

The elegant idea: reactions belong to enzyme classes based on their structural transformation.

- A **hydrolase** uses water to cleave a bond
- A **transferase** moves a group from donor to acceptor
- A **lyase** eliminates a group to form a double bond
- An **oxidoreductase** transfers electrons

This should be capturable with simple, parsimonious SMARTS patterns and stoichiometric rules. The classification should emerge naturally from the chemistry.

## What We Actually Built

In practice, the classifiers accumulate exclusions:

```python
# Transferase has 20+ specific exclusions:
if has_gtp and has_co2: return False  # "lyase"
if has_gtp and has_coa and has_pi: return False  # "ligase"
if has_nad and has_nadh: return False  # "oxidoreductase"
if single_coa_substrate: return False  # "isomerase"
# ... and so on
```

This is not parsimonious. It's a pile of patches encoding special cases.

## Why This Happens

### 1. GO/EC Classification Encodes Mechanism and Function, Not Just Structure

Consider **PEP carboxykinase**: `OAA + GTP → PEP + GDP + CO2`

Structurally, this IS a phosphoryl transfer (GTP→GDP). But it's classified as EC 4.1.1.32 (lyase) because the *biological purpose* is decarboxylation—the phosphoryl transfer is just how the enzyme achieves it mechanistically.

We're trying to match classifications that encode **intent** using only **structure**.

### 2. The Same Structural Pattern Can Mean Different Things

`ATP → ADP + Pi` appears in:

| Enzyme Type | What Happens to Phosphate | Classification |
|-------------|--------------------------|----------------|
| Kinase | Transferred to substrate | EC 2.7 (Transferase) |
| Transport ATPase | Released, energy drives transport | EC 7.x (Translocase) |
| Ligase | Released, energy forms new bond | EC 6.x (Ligase) |
| ATPase | Just hydrolyzed | EC 3.6 (Hydrolase) |

**However**, some of these ARE distinguishable with available data:

- **Transport reactions** have location annotations in RHEA: `sulfate(out) + ATP → sulfate(in) + ADP + Pi`
- **Ligases** often show `ATP → AMP + PPi` (different products)
- **Kinases** show a phosphorylated product: `substrate + ATP → substrate-P + ADP`

The issue is that our ETL doesn't fully extract all available information. For example, RHEA labels contain `(in)`/`(out)` compartment annotations, but these aren't being parsed into the `location` field of Participant objects. This is fixable.

**What truly requires external knowledge:**
- Distinguishing a kinase from an ATPase when both show `ATP → ADP + Pi` and neither has location or phosphorylated-product annotations
- Protein family, active site architecture, metabolic context

### 3. Some Ground Truth Is Probably Wrong

Examining false positives and negatives reveals suspicious annotations:

- RHEA:15040 (SAM methyltransferase) mapped to GO:0016798 (glycosidase)?
- RHEA:17368 (Mn2+ transport ATPase) mapped to GO:0016798 (glycosidase)?

We may be overfitting to incorrect GO→RHEA mappings.

### 4. Classification Is Historically Contingent

EC numbers were assigned over decades by different committees. Some boundaries are arbitrary or reflect the order of discovery rather than fundamental biochemistry. The "lyase" vs "transferase" distinction for PEP carboxykinase is a judgment call that could reasonably have gone either way.

## What Knowledge Are Enzymologists Using?

When an enzymologist classifies a reaction, they consider:

| Information Type | Example | Do We Have It? |
|-----------------|---------|----------------|
| Protein fold/family | "Rossmann fold → likely dehydrogenase" | No |
| Active site residues | "Catalytic triad → serine protease" | No |
| Cofactor chemistry | "FAD + O2 → oxidase; FAD alone → dehydrogenase" | Partial |
| Metal requirements | "Zn2+ in active site → metalloprotease" | No |
| Metabolic context | "In gluconeogenesis pathway" | No |
| Kinetic mechanism | "Ping-pong vs sequential" | No |

We have only substrate/product structures and cofactor presence. We're trying to infer mechanism from correlation, not causation.

## The Fundamental Tension

Enzyme classification conflates three distinct questions:

1. **What** chemical transformation occurs? (structure)
2. **How** does the enzyme make it happen? (mechanism)
3. **Why** does the enzyme exist? (biological function)

EC/GO classifications mix all three, but lean heavily toward mechanism and function. Our structural approach only captures #1.

## Paths Forward

### Option 1: Accept the Heuristic Nature

Acknowledge that we're building a useful approximation, not a principled classification. The rules encode "what usually correlates with X" rather than "what defines X." This is honest and may be sufficient for many applications.

**Pros**: Pragmatic, achievable
**Cons**: Won't match expert classification reliably

### Option 2: Incorporate Mechanistic Features

Add information that enzymologists actually use:

- Protein family (from UniProt)
- Active site annotations
- Cofactor binding modes (not just presence)
- Known mechanism patterns

**Pros**: More principled, higher ceiling
**Cons**: Requires additional data sources, more complex

### Option 3: Question the Target

Perhaps GO/EC isn't the right ground truth for structure-based classification. Alternatives:

- Predict EC number from reaction structure (EC is more mechanistic than GO)
- Define our own structure-based ontology
- Focus on specific sub-problems where structure IS determinative

**Pros**: Intellectually honest alignment of features and targets
**Cons**: May not match existing enzyme databases

### Option 4: Hybrid Approach

Use structural rules as a first pass, then refine with mechanism-aware rules for ambiguous cases. Accept that some reactions require additional context to classify correctly.

## Conclusions

The situation is more nuanced than "structure alone can't work":

### What IS achievable with better data extraction

Many apparent "same structure, different class" problems are solvable:

1. **Transport vs kinase**: Use location annotations `(in)`/`(out)` already in RHEA
2. **Ligase vs kinase**: Check for AMP+PPi vs ADP products
3. **Transferase vs lyase with same cofactor**: Check if a group actually transfers to an acceptor

These require **better ETL**, not external knowledge. Current gaps:
- Location annotations not parsed from RHEA labels
- Incomplete ChEBI ID coverage (many participants have null ChEBI)
- Polymer reactions lose structural detail

### What truly requires external knowledge

Some distinctions cannot be made from reaction structure alone:

1. **Same products, different mechanism**: `ATP → ADP + Pi` where phosphate isn't visibly attached to anything (true ATPase vs kinase with undetected phospho-product)
2. **Classification by biological role**: PEP carboxykinase classified as lyase despite phosphoryl transfer, because decarboxylation is the "point"
3. **Evolutionarily-based groupings**: Enzyme families grouped by ancestry, not chemistry

### Practical implications

The accumulating exclusion rules fall into two categories:

1. **Compensating for ETL gaps** (fixable): "Exclude transport if same molecule on both sides" - should use location field instead
2. **Encoding mechanistic knowledge** (inherent): "GTP + CO2 = lyase even though GTP→GDP" - genuine ambiguity

A realistic goal: **maximize what's achievable from reaction structure**, while being clear about the remaining cases that require protein-level information.

Current ~80% accuracy with interpretable rules has value. The path forward is:
1. Fix ETL to extract all available RHEA data (locations, complete ChEBI mapping)
2. Accept that some residual cases need mechanism/family information
3. Document which rules are principled vs compensatory

## Appendix: Known ETL Gaps

### Location annotations not extracted

RHEA labels contain compartment information:
```
"sulfate(out) + ATP + H2O = sulfate(in) + ADP + phosphate"
```

The `parse_location_from_label()` function exists in `rhea_etl.py`, but:
- `ParsedParticipant` class lacks a `location` field
- `parsed_to_participant()` doesn't set `Participant.location`
- Result: all `location` fields are null in cached data

**Fix**: Add location to ParsedParticipant, wire up the parsing.

### Missing ChEBI IDs

Many RHEA participants have `chebi_id: null` even for common molecules like GTP/GDP. This forces SMILES-based detection which is error-prone.

**Potential fix**: Cross-reference with ChEBI database using SMILES/InChI matching.

### Polymer notation loses structure

Reactions like `(1->4)-alpha-D-glucan(n) + H2O = (1->4)-alpha-D-glucan(n-1) + glucose` have useful stoichiometric information that could inform classification.

## References

- [EC Number Classification](https://www.qmul.ac.uk/sbcs/iubmb/enzyme/)
- [Gene Ontology Molecular Function](http://geneontology.org/docs/ontology-documentation/)
- [Rhea Reaction Database](https://www.rhea-db.org/)
