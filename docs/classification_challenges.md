# Classification Challenges v2: Strategic Direction

*Updated after reaching 77 classifiers with 1841 TP, 310 FP, micro-F1 0.788*

## Executive Summary

Autarch has evolved from a proof-of-concept into a production-grade enzyme reaction classification system. After extensive development, we've reached a decision point: **the current approach has clear strengths and well-defined limits**. This document proposes a strategic pivot that leverages what works while acknowledging what doesn't.

## Where We Are: Current Performance Profile

### By-the-Numbers

| Metric | Value |
|--------|-------|
| Total classifiers | 77 |
| Perfect classifiers (F1=1.0) | 5 |
| Good classifiers (F1≥0.7) | 20 |
| Poor classifiers (F1>0) | 23 |
| Zero-performance classifiers | 29 |
| Micro-precision | 0.856 |
| Micro-recall | 0.729 |
| Micro-F1 | 0.788 |

### What Works Well

**Perfect or near-perfect classifiers (F1 ≥ 0.9):**
- tRNASynthetase (1.0), Sulfotransferase (1.0), Peroxidase (1.0), CarbonicAnhydrase (1.0), BiotinProteinLigase (1.0)
- Isomerase (0.96), Kinase (0.93), Methyltransferase (0.92), Oxidoreductase (0.92)

These share common traits:
1. **Distinctive cofactors or products**: tRNA, biotin, peroxide, CO2
2. **Characteristic substrate classes**: kinase phosphorylates, methyltransferase uses SAM
3. **Clear stoichiometric patterns**: isomerase is 1→1

**New high-precision classifiers (recent work):**
- Dioxygenase: 8 TP, 0 FP (1.000 precision)
- CarbonSulfurLyase: 9 TP, 0 FP (1.000 precision)
- AmmoniaLyase: 5 TP, 0 FP (1.000 precision)
- CHNHOxidoreductase: 14 TP, 0 FP (1.000 precision)

### What Doesn't Work

**Zero or near-zero performance (29 classifiers with F1=0):**
- Many are highly specific (Chitinase, Xylanase, Amylase) with few or no reactions in test set
- Some are genuinely hard (Translocase requires location data, Helicase requires DNA context)
- Some have wrong GO mappings in ground truth

**Persistent challenges:**
- AldehydeOxidoreductase: 2 TP, 6 FP, 26 FN (F1=0.11)
- Racemase: 7 TP, 13 FP, 5 FN (F1=0.44)
- Reductase: 2 TP, 5 FP, 0 FN (F1=0.44)

## Root Cause Analysis

### The Data Quality Problem (Bigger Than We Thought)

1. **Sparse ground truth**: Most GO terms have only 1-4 reactions in RHEA. Hard to train or evaluate.

2. **GO term granularity mismatch**: We classify at EC level 2-3, but GO terms are often more specific (EC level 4+). A classifier for GO:0016639 (amino acid dehydrogenases) gets "false positives" that are actually children like GO:0033735 (aspartate dehydrogenase).

3. **Missing ChEBI IDs**: Many RHEA participants have null ChEBI, forcing brittle name-matching.

4. **Missing location annotations**: Translocases need `(in)`/`(out)` but RHEA label parsing isn't extracting these.

5. **Polymer/generic reactions**: "an L-alpha-amino acid" or "carbohydrate polymer" lose structural specificity.

### The Fundamental Ontology Problem

EC/GO conflate three questions:
1. **What** transformation occurs? (structure) ← We can detect this
2. **How** does the enzyme work? (mechanism) ← Partially detectable from cofactors
3. **Why** does it exist? (biological function) ← Cannot detect from structure

Example: PEP carboxykinase (`OAA + GTP → PEP + GDP + CO2`) IS a phosphotransfer structurally, but is classified as lyase (EC 4.1.1) because the "point" is decarboxylation.

## Strategic Options

### Option A: Continue Depth-First (Current Path)

Keep adding classifiers and refining existing ones.

**Pros:**
- Incremental progress is guaranteed
- Each classifier is interpretable
- Good for specific high-value targets

**Cons:**
- Diminishing returns: easy classifiers are done
- 29 zero-performance classifiers suggest limits
- Doesn't solve structural problems

**Verdict**: Continue for high-value targets, but not as primary strategy.

### Option B: Fix the Data Layer

Invest in ETL improvements:
1. Parse location annotations from RHEA labels
2. Improve ChEBI ID coverage (cross-reference with SMILES/InChI)
3. Handle polymer stoichiometry (n) → (n-1)
4. Build richer test sets from EC→RHEA mappings

**Pros:**
- Unlocks Translocase and transport-dependent classifiers
- Reduces brittle name-matching
- One-time investment with broad benefits

**Cons:**
- Significant engineering effort
- Some data simply isn't in RHEA
- Doesn't solve the ontology mismatch

**Verdict**: High priority. Location parsing alone could enable ~5-10 new classifiers.

### Option C: Embrace Hierarchical Classification

Stop treating EC level 2/3 classifiers as independent. Build explicit hierarchy:

```
Oxidoreductase (EC 1)
├── ActingOnCHOH (EC 1.1)
│   ├── WithNAD (EC 1.1.1) → AlcoholDehydrogenase
│   └── WithO2 (EC 1.1.3) → AlcoholOxidase
├── ActingOnCHNH2 (EC 1.4)
│   ├── WithNAD (EC 1.4.1) → AminoAcidDehydrogenase
│   └── WithO2 (EC 1.4.3) → MonoamineOxidase
└── ...
```

**Implementation:**
1. Parent classifies first (Oxidoreductase detects redox)
2. Children refine based on substrate/acceptor
3. Evaluation compares to appropriate GO level

**Pros:**
- Matches EC structure
- Reuses logic (parent catches general pattern, children specialize)
- Better handles GO term granularity

**Cons:**
- Requires rethinking evaluation (which level counts?)
- More complex code structure
- Some EC boundaries are arbitrary

**Verdict**: Strong conceptual fit. Would clean up much confusion.

### Option D: Pivot to Explanation, Not Classification

What if the goal isn't perfect classification, but **interpretable explanation**?

Use case: Given a reaction, explain WHY it might be classified as X.

```python
result = classifier.explain(reaction)
# Returns:
# - "Has NAD+/NADH redox pair → likely oxidoreductase"
# - "Phosphoryl group transfer from ATP → likely kinase or ligase"
# - "Uses H2O to cleave bond → consistent with hydrolase"
# - "Confidence: HIGH for oxidoreductase, MEDIUM for dehydrogenase, LOW for kinase"
```

**Pros:**
- Doesn't need perfect accuracy to be useful
- Explanations help human experts
- Naturally handles ambiguous cases
- Could be used to QC existing annotations

**Cons:**
- Requires good UX design
- "Explanation" quality is subjective
- Harder to benchmark

**Verdict**: Compelling for real-world applications. Partially implemented already (ClassificationResult.explanation).

### Option E: Machine Learning Hybrid

Train ML models on RHEA reactions, using structural features:
1. Morgan fingerprints of reactants/products
2. Change in functional groups
3. Cofactor presence vectors
4. Bond difference analysis

Use rule-based classifiers as features or for interpretability.

**Pros:**
- Can learn patterns we miss
- Handles feature interactions
- Could achieve higher recall

**Cons:**
- Less interpretable
- Requires more data than we have for many classes
- Deployment complexity

**Verdict**: Worth exploring as ensemble with rule-based. Not primary strategy.

## Recommended Path Forward

### Phase 1: Data Layer (Immediate)

1. **Parse location annotations** from RHEA labels → enables Translocase
2. **Improve ChEBI mapping** using SMILES/InChI cross-reference
3. **Build EC-based test sets** (not just GO-linked reactions)

### Phase 2: Hierarchical Refactoring (Near-term)

1. Implement parent→child classification chain
2. Refactor Oxidoreductase family (EC 1.x.x) as proof of concept
3. Update evaluation to handle hierarchy

### Phase 3: Explanation Mode (Medium-term)

1. Enhance ClassificationResult with structured explanations
2. Build "reaction analysis" tool that shows all applicable patterns
3. Focus on QC use case: "Does this RHEA→GO mapping make sense?"

### Phase 4: Selective Expansion (Ongoing)

Continue adding classifiers where:
- Ground truth has ≥10 reactions
- Structural pattern is distinctive
- Use case demands it

## Appendix: Technical Debt

### Code Issues
- Many classifiers have ad-hoc exclusions that should be pushed to parent
- Pattern DSL is underutilized (most classifiers are imperative Python)
- No test coverage for many classifiers

### Architecture Issues
- Classifier discovery is implicit (imports in `__init__.py`)
- No validation that GO_ID/EC_NUMBER_PREFIX are consistent
- Evaluation doesn't account for GO hierarchy

### Documentation Issues
- No standardized docstring format for classifiers
- Change history only in some files
- No catalog of "known hard cases"

## Conclusions

The core insight: **we've proven that structural rules work for ~40% of enzyme classes with >70% F1, and fail for the rest**. The path forward isn't more of the same—it's:

1. **Fix data quality** to unlock currently-impossible classifiers
2. **Embrace hierarchy** to handle GO granularity
3. **Pivot to explanation** for real-world value
4. **Accept limits** for mechanism-dependent classifications

The 80% of value comes from 20% of classifiers (Oxidoreductase, Transferase, Hydrolase, Lyase, Isomerase, Ligase cover most reactions). Perfect classification of long-tail enzyme families isn't achievable from structure alone—and that's okay.
