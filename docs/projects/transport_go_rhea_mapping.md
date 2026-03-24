# Project: Transport Reaction GO-RHEA Mapping Gap

## Problem Statement

RHEA contains 51 transport reactions with location annotations `(in)`/`(out)`, but only 15 have GO term mappings. This gap prevents:
1. Proper evaluation of translocase classifiers in autarch
2. Computational discovery of transport reactions via GO
3. Integration of transport knowledge across databases

## Current State (After ETL Fix)

| Metric | Before Fix | After Fix |
|--------|-----------|-----------|
| Transport reactions in RHEA (with location) | 51 | 97 |
| With GO mappings | 15 (29%) | 87 (90%) |
| Without GO mappings | 36 | 10 |
| Truly unmapped (no EC→GO) | 4 | 3 |

## Root Cause Analysis

### Original Issues (Now Fixed)

1. **ETL bug: Wrong master ID formula** - The code used `(id // 4) * 4` to compute
   RHEA master IDs, but RHEA master IDs are NOT always divisible by 4. They're
   defined in rhea-directions.tsv. This caused ~3600 GO→RHEA mappings to be lost.

2. **ETL bug: EC-only filtering** - The RHEA cache only loaded reactions with EC
   annotations, missing many reactions that have GO mappings but no EC numbers.

### Resolution

Fixed in `rhea_tsv_etl.py`:
- Use RHEA IDs directly from GO (source of truth) instead of computing master IDs
- Include reactions with GO mappings OR EC annotations (not just EC)

### Remaining Gaps

1. **EC 7.x.x.x is new** (added 2018) - some EC numbers still lack GO terms
2. **3 reactions truly unmapped** - no GO terms have their EC numbers:
   - RHEA:29791 (glutathione, EC 7.4.2.10)
   - RHEA:30259 (menaquinol:O2, EC 7.1.1.5)
   - RHEA:47748 (NADH:ubiquinone, EC 7.2.1.1)

## Proposed Solution

### Phase 1: Direct EC-based Mapping

Many unmapped RHEA reactions have EC 7.x.x.x numbers. GO terms often include EC numbers.

**Strategy**: Match RHEA reactions to GO terms via shared EC numbers.

Example:
- RHEA:14633 has EC 7.2.2.3, 7.2.2.4 (sodium ATPases)
- GO:0015081 "sodium ion transmembrane transporter activity" should link

### Phase 2: Substrate-based Mapping

For reactions without EC matches, use transported substrate:
- Mg2+ transport → GO terms containing "magnesium"
- Ca2+ transport → GO terms containing "calcium"

### Phase 3: Manual Curation

Review generated mappings and flag uncertain cases for expert review.

## Deliverables

1. `proposed_mappings.tsv` - Machine-generated RHEA→GO proposals
2. `mapping_evidence.tsv` - Evidence supporting each mapping
3. `uncertain_mappings.tsv` - Cases needing manual review
4. Summary statistics

## Implementation

See: `scripts/generate_transport_mappings.py`

## Timeline

- Phase 1: Automated generation (this session)
- Phase 2: Review and refinement
- Phase 3: Submission to RHEA/GO teams

## Success Criteria

- Generate high-confidence mappings for ≥25 of 36 unmapped reactions
- <5% error rate on spot-checked mappings
- Mappings accepted by RHEA/GO curators

---

## Execution Log

### 2024-12-27: Initial Analysis

**Unmapped transport reactions (36 total):**

| RHEA ID | Transported | EC Numbers | Status |
|---------|-------------|------------|--------|
| RHEA:13181 | nitrate | 7.3.2.4 | Needs mapping |
| RHEA:13973 | H+ | 7.1.3.1 | Needs mapping |
| RHEA:14613 | taurine | 7.6.2.7 | Needs mapping |
| RHEA:14633 | Na+ | 7.2.2.3, 7.2.2.4 | Needs mapping |
| RHEA:14733 | Ag+ | 7.2.2.15 | Needs mapping |
| RHEA:15557 | Ni2+ | 7.2.2.11 | Needs mapping |
| RHEA:16777 | K+ | 7.2.2.6 | Needs mapping |
| RHEA:17365 | Mn2+ | 7.2.2.5 | Needs mapping |
| RHEA:18065 | phosphonate | 7.3.2.2 | Needs mapping |
| RHEA:18105 | Ca2+ | 7.2.2.10 | Needs mapping |
| RHEA:18353 | Na+/K+ | 7.2.2.13 | Needs mapping |
| RHEA:19261 | heme b | 7.6.2.5 | Needs mapping |
| RHEA:20621 | Zn2+ | 7.2.2.12 | Needs mapping |
| RHEA:21396 | Na+ | 7.2.4.3 | Needs mapping |
| RHEA:29767 | D-methionine | 7.4.2.11 | Needs mapping |
| RHEA:29779 | L-methionine | 7.4.2.11 | Needs mapping |
| RHEA:29791 | glutathione | 7.4.2.10 | Needs mapping |
| RHEA:29795 | Zn2+ | 7.2.2.20 | Needs mapping |
| RHEA:29799 | D-allose | 7.5.2.8 | Needs mapping |
| RHEA:29811 | thiamine | 7.6.2.15 | Needs mapping |
| RHEA:29871 | thiosulfate | 7.3.2.3 | Needs mapping |
| RHEA:29879 | L-arginine | 7.4.2.1 | Needs mapping |
| RHEA:29883 | L-ornithine | 7.4.2.1 | Needs mapping |
| RHEA:29887 | L-lysine | 7.4.2.1 | Needs mapping |
| RHEA:29899 | D-xylose | 7.5.2.10, 7.5.2.13 | Needs mapping |
| RHEA:29903 | D-ribose | 7.5.2.7 | Needs mapping |
| RHEA:29995 | putrescine | 7.6.2.16 | Needs mapping |
| RHEA:29999 | spermidine | 7.6.2.11 | Needs mapping |
| RHEA:30007 | L-arabinose | 7.5.2.12, 7.5.2.13 | Needs mapping |
| RHEA:33147 | daunorubicin | 7.6.2.2 | Needs mapping |
| RHEA:35027 | tungstate | 7.3.2.6 | Needs mapping |
| RHEA:36831 | H+ | 7.1.1.4 | Needs mapping |
| RHEA:30259 | H+ | 7.1.1.5 | Needs mapping |
| RHEA:47748 | Na+ | 7.2.1.1 | Needs mapping |
| RHEA:30251 | H+ | 7.1.1.3 | Needs mapping |
| RHEA:43336 | Na+ | 7.2.4.1 | Needs mapping |

**GO terms that should have RHEA links but don't:**
- GO:0015081 (sodium ion transmembrane transporter)
- GO:0015085 (calcium ion transmembrane transporter)
- GO:0015079 (potassium ion transmembrane transporter)
- GO:0015099 (nickel cation transmembrane transporter)
- GO:0005385 (zinc ion transmembrane transporter)
- GO:0015095 (magnesium ion transmembrane transporter)
- ... and many more

---

## Execution Results (2024-12-27)

### ETL Fix Impact

The ETL fix dramatically improved GO coverage:

| Metric | Before Fix | After Fix |
|--------|-----------|-----------|
| Total RHEA reactions cached | ~4760 | 4925 |
| Reactions with GO mappings | 1389 (28%) | 3814 (77%) |
| Transport reactions | 51 | 97 |
| Transport with GO mappings | 15 (29%) | 87 (90%) |

### Remaining Proposed Mappings

After the fix, only 10 transport reactions lack GO mappings:

| Category | Count |
|----------|-------|
| High confidence (EC match) | 6 |
| Medium confidence (substrate match) | 54 |
| RHEA reactions covered | 7 / 10 |
| RHEA reactions still unmapped | 3 |

### Still Unmapped (3 reactions)

These reactions lack matching GO terms with EC numbers:

1. **RHEA:29791** - glutathione transport (EC 7.4.2.10)
   - No GO term has this EC number

2. **RHEA:30259** - menaquinol:O2 oxidoreductase with H+ transport (EC 7.1.1.5)
   - Complex electron transport chain reaction

3. **RHEA:47748** - NADH:ubiquinone with Na+ transport (EC 7.2.1.1)
   - No GO term has this EC number

### Files Generated

1. `proposed_mappings.tsv` - All mappings with confidence levels
2. `rhea_go_proposed_mappings_for_submission.tsv` - Clean submission file (high-confidence only)

### Next Steps

1. **Review** - Manual review of high-confidence mappings
2. **Submit to RHEA** - Contact RHEA curators with proposed mappings
3. **Submit to GO** - Or contact GO curators to add RHEA cross-references
4. **Follow up on missing EC→GO** - The 4 unmapped reactions indicate missing GO terms

### Impact

If these mappings are accepted:
- 31 additional RHEA reactions will have GO annotations
- autarch Translocase classifier can be re-enabled
- Transport reaction classification accuracy will improve significantly
