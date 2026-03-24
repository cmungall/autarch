# ETL Improvement Plan

Based on exhaustive analysis of the RHEA data pipeline, this document outlines specific gaps and fixes.

## Executive Summary

| Issue | Impact | Effort | Priority |
|-------|--------|--------|----------|
| Location annotations not extracted | 42 transport reactions unidentifiable | Low | **High** |
| ChEBI ID normalization | 263 common molecules missing IDs | Low | **High** |
| RHEA EQUATION parsing | Could fill more ChEBI IDs | Medium | Medium |
| Polymer stoichiometry handling | 218 participants affected | Medium | Low |

## Detailed Analysis

### 1. Location Annotations (HIGH PRIORITY)

**Current State:**
- 42 reactions have `(in)`/`(out)` annotations in RHEA labels
- 56 occurrences each of `(in)` and `(out)`
- `Participant.location` field exists but is always `null`
- `parse_location_from_label()` function exists in `rhea_etl.py`

**Example:**
```
Label: "sulfate(out) + ATP + H2O = sulfate(in) + ADP + phosphate + H(+)"
Current: location=null for all participants
Should be: sulfate has location="out" on left, location="in" on right
```

**Root Cause:**
- `ParsedParticipant` class lacks `location` field
- `parsed_to_participant()` doesn't set location
- Location parsing exists but isn't wired up

**Fix:**
1. Add `location: Optional[str] = None` to `ParsedParticipant` in `label_parser.py`
2. Update label parser to extract `(in)`, `(out)`, `(cytoplasm)`, etc.
3. Update `parsed_to_participant()` to set `Participant.location`
4. Regenerate cache

**Files to modify:**
- `src/autarch/etl/label_parser.py`
- `src/autarch/etl/rhea_tsv_etl.py`

---

### 2. ChEBI ID Normalization (HIGH PRIORITY)

**Current State:**
- 2,570 participants (18.8%) lack ChEBI IDs
- 263 of these are common molecules that SHOULD have ChEBI IDs
- Same molecule appears with ChEBI in some reactions, without in others

**Molecules affected:**
| Molecule | Missing Count | Correct ChEBI |
|----------|---------------|---------------|
| phosphate | 83 | CHEBI:43474 |
| H(+) | 68 | CHEBI:15378 |
| H2O | 23 | CHEBI:15377 |
| UDP | 18 | CHEBI:58223 |
| sulfur | 15 | CHEBI:17909 |
| O2 | 14 | CHEBI:15379 |
| ADP | 10 | CHEBI:456216 |
| hydrogen sulfide | 9 | CHEBI:16136 |
| ATP | 8 | CHEBI:30616 |

**Root Cause:**
- Polymer reactions with `n` stoichiometry lose ChEBI IDs
- Some RHEA entries inconsistently annotated

**Fix:**
1. Create normalization dictionary mapping names → ChEBI IDs
2. Apply during ETL when ChEBI ID is null but name matches
3. Consider SMILES-based matching as fallback

**Implementation:**
```python
CHEBI_BY_NAME = {
    'H(+)': 'CHEBI:15378',
    'hydron': 'CHEBI:15378',
    'water': 'CHEBI:15377',
    'H2O': 'CHEBI:15377',
    'phosphate': 'CHEBI:43474',
    # ... etc
}

def normalize_chebi(participant):
    if participant.chebi_id:
        return participant
    name = participant.name.lower().strip()
    if name in CHEBI_BY_NAME:
        participant.chebi_id = CHEBI_BY_NAME[name]
    return participant
```

---

### 3. RHEA EQUATION Field (MEDIUM PRIORITY)

**Current State:**
- `rhea-reactions.txt.gz` contains EQUATION field with explicit ChEBI IDs
- Example: `CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938`
- This is not currently used

**Potential:**
- Could cross-reference to fill ChEBI IDs where SMILES matching fails
- More reliable than name matching

**Fix:**
1. Parse EQUATION field from RHEA text file
2. Map ChEBI IDs to participant positions
3. Use as authoritative source for ChEBI assignment

---

### 4. Polymer Stoichiometry (LOW PRIORITY)

**Current State:**
- 218 participants (1.6%) have variable stoichiometry (`n`, `2n`, `n+1`, `n-1`)
- These often lack ChEBI IDs and SMILES
- Example: `n malonyl-CoA + acetyl-CoA + 2n NADPH + 2n H(+)`

**Impact:**
- Affects classification accuracy for polymer reactions
- Loss of quantitative information

**Fix:**
- Mostly a display/representation issue
- Could normalize common molecules (H+, NADPH) even in polymer context
- Consider flagging polymer reactions for special handling

---

## Breakdown of Missing ChEBI IDs

Total missing: **2,570 participants**

| Category | Count | % | Recoverable? |
|----------|-------|---|--------------|
| Generic polymer | 907 | 35% | No - by design |
| Compound by SMILES | 905 | 35% | Partial - InChI lookup |
| Protein references | 248 | 10% | No - generic |
| Common molecules | 263 | 10% | **Yes - name normalization** |
| Other | 247 | 10% | Partial |

**Realistically recoverable: ~300 participants (12%)**

---

## Implementation Order

### Phase 1: Quick Wins (1-2 hours)
1. Add ChEBI normalization by name (263 fixes)
2. Regenerate cache
3. Re-run evaluation

### Phase 2: Location Annotations (2-4 hours)
1. Add location field to ParsedParticipant
2. Wire up location parsing
3. Regenerate cache
4. Update classifiers to use location info

### Phase 3: Advanced (Optional)
1. RHEA EQUATION parsing for ChEBI recovery
2. InChI-based ChEBI lookup
3. Polymer reaction flagging

---

## Impact on Classification

With these fixes:

| Classifier | Current Issue | After Fix |
|------------|---------------|-----------|
| Transferase | Transport detection via heuristics | Use `location` field |
| Kinase | ATP→ADP heuristics | Use `location` for transport exclusion |
| ATPase | Overlap with kinases | Location distinguishes transport |
| All | Missing cofactor detection | Better ChEBI coverage |

**Estimated improvement:**
- Transport vs kinase confusion: **fully resolvable**
- Missing cofactor patterns: **10% more detectable**
- Overall: **2-5% F1 improvement possible**

---

## Validation

After implementing fixes:
1. Count participants with location annotations (expect 56+)
2. Count participants with ChEBI IDs (expect 81.2% → ~83%)
3. Re-run `just eval-all` and compare metrics
4. Specifically check transport reaction classification
