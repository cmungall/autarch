#!/usr/bin/env python3
"""Generate proposed GO-RHEA mappings for transport reactions.

This script identifies unmapped transport reactions in RHEA and proposes
GO term mappings based on:
1. Shared EC numbers
2. Substrate/transported molecule matching
3. GO term label matching

When a GO term has a mechanistic qualifier (ABC-type, P-type, etc.),
the script looks for parent GO terms without the qualifier to provide
more accurate exactMatch mappings.
"""

import json
import re
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class ProposedMapping:
    rhea_id: str
    rhea_label: str
    go_id: str
    go_label: str
    predicate_id: str  # skos:exactMatch, skos:broadMatch, etc.
    evidence_type: str  # "ec_match", "substrate_match", "label_match"
    evidence_detail: str
    confidence: str  # "high", "medium", "low"


# Mechanistic qualifiers that make GO terms more specific than RHEA
MECHANISTIC_QUALIFIERS = [
    'abc-type',
    'p-type',
    'v-type',
    'f-type',
    'secondary active',
    'primary active',
    'symporter',
    'antiporter',
    'uniporter',
]


def load_data(cache_dir: str = "cache"):
    """Load RHEA reactions and GO terms."""
    rhea_reactions = []
    with open(Path(cache_dir) / "rhea_reactions.jsonl") as f:
        for line in f:
            rhea_reactions.append(json.loads(line))

    go_terms = {}
    with open(Path(cache_dir) / "go_terms.jsonl") as f:
        for line in f:
            t = json.loads(line)
            go_terms[t["go_id"]] = t

    return rhea_reactions, go_terms


def has_mechanistic_qualifier(label: str) -> tuple[bool, str | None]:
    """Check if GO label has a mechanistic qualifier."""
    label_lower = label.lower()
    for qualifier in MECHANISTIC_QUALIFIERS:
        if qualifier in label_lower:
            return True, qualifier
    return False, None


def find_parent_without_mechanism(go_id: str, go_terms: dict, substrate_keywords: list[str]) -> tuple[str, str] | None:
    """Find an ancestor GO term that matches substrate but lacks mechanism qualifier.

    Returns (go_id, label) of best parent, or None if not found.
    """
    go_term = go_terms.get(go_id)
    if not go_term:
        return None

    ancestors = go_term.get("ancestors", [])
    candidates = []

    for ancestor_id in ancestors:
        if ancestor_id == go_id:  # Skip self
            continue
        ancestor = go_terms.get(ancestor_id)
        if not ancestor:
            continue

        ancestor_label = ancestor.get("label", "").lower()

        # Check if ancestor has substrate match
        has_substrate = any(kw in ancestor_label for kw in substrate_keywords)
        if not has_substrate:
            continue

        # Check if ancestor lacks mechanistic qualifier
        has_mech, _ = has_mechanistic_qualifier(ancestor_label)
        if has_mech:
            continue

        # Check it's a transporter term
        if "transport" not in ancestor_label:
            continue

        # Score by specificity (longer label = more specific, prefer that)
        candidates.append((ancestor_id, ancestor.get("label", ""), len(ancestor_label)))

    if not candidates:
        return None

    # Return most specific (longest label) parent without mechanism
    candidates.sort(key=lambda x: -x[2])
    return (candidates[0][0], candidates[0][1])


def extract_substrate_keywords(rhea_label: str) -> list[str]:
    """Extract substrate keywords from RHEA reaction label."""
    keywords = []

    # Ion patterns
    ion_patterns = {
        r'\bna\b|\bsodium\b|\bna\(': ['sodium', 'na'],
        r'\bk\b|\bpotassium\b|\bk\(': ['potassium', 'k+'],
        r'\bca\b|\bcalcium\b|\bca\(': ['calcium', 'ca'],
        r'\bmg\b|\bmagnesium\b|\bmg\(': ['magnesium', 'mg'],
        r'\bzn\b|\bzinc\b|\bzn\(': ['zinc', 'zn'],
        r'\bfe\b|\biron\b|\bfe\(': ['iron', 'fe'],
        r'\bcu\b|\bcopper\b|\bcu\(': ['copper', 'cu'],
        r'\bmn\b|\bmanganese\b|\bmn\(': ['manganese', 'mn'],
        r'\bni\b|\bnickel\b|\bni\(': ['nickel', 'ni'],
        r'\bag\b|\bsilver\b|\bag\(': ['silver', 'ag'],
        r'\bh\b|\bproton\b|\bhydron\b|\bh\(': ['proton', 'h+'],
        r'\bnitrate\b': ['nitrate'],
        r'\bsulfate\b|\bsulphate\b': ['sulfate'],
        r'\bphosphate\b': ['phosphate'],
        r'\bphosphonate\b': ['phosphonate'],
        r'\bthiosulfate\b': ['thiosulfate'],
        r'\btungstate\b': ['tungstate'],
    }

    # Organic molecule patterns
    organic_patterns = {
        r'\btaurine\b': ['taurine'],
        r'\bmethionine\b': ['methionine'],
        r'\barginine\b': ['arginine'],
        r'\blysine\b': ['lysine'],
        r'\bornithine\b': ['ornithine'],
        r'\bglutathione\b': ['glutathione'],
        r'\bxylose\b': ['xylose'],
        r'\bribose\b': ['ribose'],
        r'\barabinose\b': ['arabinose'],
        r'\ballose\b': ['allose'],
        r'\bthiamine\b': ['thiamine'],
        r'\bputrescine\b': ['putrescine'],
        r'\bspermidine\b': ['spermidine', 'polyamine'],
        r'\bheme\b': ['heme'],
        r'\bdaunorubicin\b': ['xenobiotic', 'drug'],
    }

    label_lower = rhea_label.lower()

    for pattern, kws in {**ion_patterns, **organic_patterns}.items():
        if re.search(pattern, label_lower):
            keywords.extend(kws)

    return keywords


def determine_predicate_and_target(
    rhea_id: str,
    rhea_label: str,
    go_id: str,
    go_label: str,
    go_terms: dict,
    ec_numbers: list[str]
) -> tuple[str, str, str, str]:
    """Determine SKOS predicate and potentially find better GO target.

    Returns (target_go_id, target_go_label, predicate_id, justification).
    """
    go_lower = go_label.lower()
    rhea_lower = rhea_label.lower()

    # Check for mechanistic qualifier
    has_mech, qualifier = has_mechanistic_qualifier(go_label)

    if has_mech:
        # Try to find parent without mechanism
        substrate_keywords = extract_substrate_keywords(rhea_label)
        parent = find_parent_without_mechanism(go_id, go_terms, substrate_keywords)

        if parent:
            parent_id, parent_label = parent
            ec_str = ', '.join(ec_numbers) if ec_numbers else 'unknown'
            return (
                parent_id,
                parent_label,
                'skos:exactMatch',
                f"Shared EC: {ec_str}; mapped_to_parent_without_mechanism"
            )
        else:
            # No suitable parent found - use broadMatch to mechanism-specific term
            qualifier_clean = qualifier.replace('-', '_').replace(' ', '_')
            ec_str = ', '.join(ec_numbers) if ec_numbers else 'unknown'
            return (
                go_id,
                go_label,
                'skos:broadMatch',
                f"Shared EC: {ec_str}; go_has_mechanism_{qualifier_clean}"
            )

    # No mechanistic qualifier - check for stereoisomer issues
    stereoisomer_keywords = ['l-', 'd-', '(s)-', '(r)-']
    rhea_has_stereo = any(s in rhea_lower for s in stereoisomer_keywords)
    go_has_stereo = any(s in go_lower for s in stereoisomer_keywords)

    if rhea_has_stereo != go_has_stereo:
        ec_str = ', '.join(ec_numbers) if ec_numbers else 'unknown'
        return (go_id, go_label, 'skos:closeMatch', f"Shared EC: {ec_str}; stereoisomer_mismatch")

    # Clean match - exactMatch
    ec_str = ', '.join(ec_numbers) if ec_numbers else 'unknown'
    return (go_id, go_label, 'skos:exactMatch', f"Shared EC: {ec_str}; substrate_match")


def get_transport_reactions(reactions):
    """Find reactions with location annotations."""
    transport = []
    for r in reactions:
        rxn = r.get("reaction", {})
        locations = []
        for side in ["left_participants", "right_participants"]:
            for p in rxn.get(side, []):
                loc = p.get("location")
                if loc:
                    locations.append({
                        "name": p.get("name", ""),
                        "location": loc,
                        "chebi_id": p.get("chebi_id")
                    })
        if locations:
            transport.append({
                "rhea_id": r["rhea_id"],
                "label": r.get("label", ""),
                "go_terms": r.get("go_terms", []),
                "ec_numbers": r.get("ec_numbers", []),
                "transported": locations
            })
    return transport


def get_transporter_go_terms(go_terms):
    """Find GO terms related to transport."""
    transporter_terms = []
    for go_id, t in go_terms.items():
        label = t.get("label", "").lower()
        if "transport" in label or "translocase" in label or "atpase" in label:
            transporter_terms.append({
                "go_id": go_id,
                "label": t.get("label", ""),
                "ec_numbers": t.get("ec_numbers", []),
                "rhea_ids": t.get("rhea_ids", []),
                "ancestors": t.get("ancestors", [])
            })
    return transporter_terms


def match_by_ec(reaction, go_terms_list, go_terms_dict):
    """Match reaction to GO terms by EC number."""
    matches = []
    rxn_ecs = set(reaction.get("ec_numbers", []))

    for go_term in go_terms_list:
        go_ecs = set(go_term.get("ec_numbers", []))
        shared = rxn_ecs & go_ecs
        if shared:
            # Determine predicate and potentially find better target
            target_id, target_label, predicate, justification = determine_predicate_and_target(
                reaction["rhea_id"],
                reaction["label"],
                go_term["go_id"],
                go_term["label"],
                go_terms_dict,
                list(shared)
            )

            matches.append(ProposedMapping(
                rhea_id=reaction["rhea_id"],
                rhea_label=reaction["label"],
                go_id=target_id,
                go_label=target_label,
                predicate_id=predicate,
                evidence_type="ec_match",
                evidence_detail=justification,
                confidence="high"
            ))
    return matches


def normalize_substrate(name):
    """Normalize substrate name for matching."""
    name = name.lower()
    # Remove charge notation
    name = re.sub(r'\(\d*[+-]\)', '', name)
    # Remove common suffixes
    name = re.sub(r' zwitterion$', '', name)
    name = re.sub(r' ion$', '', name)
    # Common normalizations
    name = name.replace("(1+)", "").replace("(2+)", "").replace("(3+)", "")
    name = name.strip()
    return name


# Substrate to GO term keyword mapping
SUBSTRATE_GO_KEYWORDS = {
    "sodium": ["sodium", "na+", "na(+)"],
    "potassium": ["potassium", "k+", "k(+)"],
    "calcium": ["calcium", "ca2+", "ca(2+)"],
    "magnesium": ["magnesium", "mg2+", "mg(2+)"],
    "zinc": ["zinc", "zn2+", "zn(2+)"],
    "iron": ["iron", "fe", "ferric", "ferrous"],
    "copper": ["copper", "cu+", "cu2+", "cu(+)", "cu(2+)"],
    "manganese": ["manganese", "mn2+", "mn(2+)"],
    "nickel": ["nickel", "ni2+", "ni(2+)"],
    "cobalt": ["cobalt", "co2+", "co(2+)"],
    "silver": ["silver", "ag+", "ag(+)"],
    "cadmium": ["cadmium", "cd2+", "cd(2+)"],
    "proton": ["proton", "h+", "h(+)", "hydron"],
    "phosphate": ["phosphate", "phosphonate"],
    "sulfate": ["sulfate", "sulphate"],
    "nitrate": ["nitrate"],
    "chloride": ["chloride", "cl-"],
    "amino acid": ["amino acid", "arginine", "lysine", "ornithine", "methionine",
                   "glutamate", "aspartate", "glycine", "alanine"],
    "sugar": ["glucose", "fructose", "maltose", "xylose", "ribose", "arabinose",
              "allose", "galactose"],
    "polyamine": ["putrescine", "spermidine", "spermine"],
}


def match_by_substrate(reaction, go_terms_list, go_terms_dict):
    """Match reaction to GO terms by transported substrate."""
    matches = []

    # Get transported molecules
    transported = [normalize_substrate(t["name"]) for t in reaction.get("transported", [])]
    transported_text = " ".join(transported)

    for go_term in go_terms_list:
        go_label = go_term["label"].lower()

        # Check each substrate category
        for category, keywords in SUBSTRATE_GO_KEYWORDS.items():
            # Does reaction transport this category?
            rxn_has = any(kw in transported_text for kw in keywords)
            # Does GO term mention this category?
            go_has = category in go_label or any(kw in go_label for kw in keywords)

            if rxn_has and go_has:
                # Determine predicate
                target_id, target_label, predicate, justification = determine_predicate_and_target(
                    reaction["rhea_id"],
                    reaction["label"],
                    go_term["go_id"],
                    go_term["label"],
                    go_terms_dict,
                    reaction.get("ec_numbers", [])
                )

                matches.append(ProposedMapping(
                    rhea_id=reaction["rhea_id"],
                    rhea_label=reaction["label"],
                    go_id=target_id,
                    go_label=target_label,
                    predicate_id=predicate,
                    evidence_type="substrate_match",
                    evidence_detail=justification,
                    confidence="medium"
                ))
                break  # One match per GO term

    return matches


def deduplicate_mappings(mappings):
    """Keep best mapping for each RHEA-GO pair."""
    best = {}
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    predicate_order = {"skos:exactMatch": 0, "skos:closeMatch": 1, "skos:broadMatch": 2, "skos:narrowMatch": 3}

    for m in mappings:
        key = (m.rhea_id, m.go_id)
        if key not in best:
            best[key] = m
        else:
            # Prefer higher confidence, then better predicate
            existing = best[key]
            if confidence_order.get(m.confidence, 99) < confidence_order.get(existing.confidence, 99):
                best[key] = m
            elif confidence_order.get(m.confidence, 99) == confidence_order.get(existing.confidence, 99):
                if predicate_order.get(m.predicate_id, 99) < predicate_order.get(existing.predicate_id, 99):
                    best[key] = m

    return list(best.values())


def filter_specific_go_terms(mappings, go_terms):
    """Prefer more specific GO terms over generic ones."""
    # Generic terms to deprioritize
    generic_terms = {
        "GO:0015075",  # monoatomic ion transmembrane transporter activity
        "GO:0022857",  # transmembrane transporter activity
        "GO:0005215",  # transporter activity
        "GO:0022890",  # inorganic cation transmembrane transporter
    }

    # Group by RHEA ID
    by_rhea = defaultdict(list)
    for m in mappings:
        by_rhea[m.rhea_id].append(m)

    filtered = []
    for rhea_id, group in by_rhea.items():
        # Prefer specific over generic
        specific = [m for m in group if m.go_id not in generic_terms]
        if specific:
            filtered.extend(specific)
        else:
            filtered.extend(group)

    return filtered


def main():
    print("Loading data...")
    reactions, go_terms = load_data()

    print("Finding transport reactions...")
    transport = get_transport_reactions(reactions)
    unmapped = [r for r in transport if not r["go_terms"]]
    print(f"  Total transport reactions: {len(transport)}")
    print(f"  Unmapped: {len(unmapped)}")

    print("Finding transporter GO terms...")
    transporter_go = get_transporter_go_terms(go_terms)
    print(f"  Found {len(transporter_go)} transporter GO terms")

    print("Generating mappings...")
    all_mappings = []

    for rxn in unmapped:
        # Try EC matching first (highest confidence)
        ec_matches = match_by_ec(rxn, transporter_go, go_terms)
        all_mappings.extend(ec_matches)

        # Try substrate matching
        substrate_matches = match_by_substrate(rxn, transporter_go, go_terms)
        all_mappings.extend(substrate_matches)

    print(f"  Raw mappings: {len(all_mappings)}")

    # Deduplicate
    all_mappings = deduplicate_mappings(all_mappings)
    print(f"  After dedup: {len(all_mappings)}")

    # Filter to prefer specific terms
    all_mappings = filter_specific_go_terms(all_mappings, go_terms)
    print(f"  After filtering: {len(all_mappings)}")

    # Sort by confidence then RHEA ID
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    all_mappings.sort(key=lambda m: (confidence_order[m.confidence], m.rhea_id))

    # Output
    output_dir = Path("docs/projects")

    # SSSOM-style TSV output
    with open(output_dir / "rhea_go_proposed_mappings_for_submission.tsv", "w") as f:
        f.write("# Proposed RHEA-GO mappings for transport reactions (SSSOM-style)\n")
        f.write("# Generated by autarch project analysis\n")
        f.write("# Predicates:\n")
        f.write("#   skos:exactMatch - substrate matches, no mechanistic difference\n")
        f.write("#   skos:broadMatch - RHEA is broader (GO has mechanism like ABC-type)\n")
        f.write("#   skos:narrowMatch - RHEA is narrower (GO is general category)\n")
        f.write("#   skos:closeMatch - related but with differences\n")
        f.write("# Date: 2024-12-27\n")
        f.write("#\n")
        f.write("subject_id\tsubject_label\tpredicate_id\tobject_id\tobject_label\tmapping_justification\n")
        for m in all_mappings:
            if m.confidence == "high":  # Only high-confidence for submission
                rhea_label_short = m.rhea_label[:80] + "..." if len(m.rhea_label) > 80 else m.rhea_label
                f.write(f"{m.rhea_id}\t{rhea_label_short}\t{m.predicate_id}\t{m.go_id}\t{m.go_label}\t{m.evidence_detail}\n")

    # Also write full proposed_mappings.tsv with all confidence levels
    with open(output_dir / "proposed_mappings.tsv", "w") as f:
        f.write("rhea_id\trhea_label\tpredicate_id\tgo_id\tgo_label\tevidence_type\tevidence_detail\tconfidence\n")
        for m in all_mappings:
            f.write(f"{m.rhea_id}\t{m.rhea_label}\t{m.predicate_id}\t{m.go_id}\t{m.go_label}\t{m.evidence_type}\t{m.evidence_detail}\t{m.confidence}\n")

    # Summary
    high_conf = [m for m in all_mappings if m.confidence == "high"]
    medium_conf = [m for m in all_mappings if m.confidence == "medium"]

    covered_rhea = set(m.rhea_id for m in all_mappings)
    uncovered = [r for r in unmapped if r["rhea_id"] not in covered_rhea]

    # Predicate distribution
    predicate_counts = defaultdict(int)
    for m in high_conf:
        predicate_counts[m.predicate_id] += 1

    print()
    print("=== SUMMARY ===")
    print(f"High confidence mappings: {len(high_conf)}")
    print(f"Medium confidence mappings: {len(medium_conf)}")
    print(f"RHEA reactions covered: {len(covered_rhea)} / {len(unmapped)}")
    print(f"RHEA reactions still unmapped: {len(uncovered)}")
    print()
    print("Predicate distribution (high confidence):")
    for pred, count in sorted(predicate_counts.items()):
        print(f"  {pred}: {count}")
    print()
    print(f"Output written to: {output_dir / 'rhea_go_proposed_mappings_for_submission.tsv'}")

    # Show sample exactMatch mappings
    print()
    print("=== SAMPLE exactMatch MAPPINGS ===")
    exact_matches = [m for m in high_conf if m.predicate_id == "skos:exactMatch"][:5]
    for m in exact_matches:
        print(f"{m.rhea_id} -> {m.go_id}")
        print(f"  RHEA: {m.rhea_label[:70]}...")
        print(f"  GO: {m.go_label}")
        print(f"  Predicate: {m.predicate_id}")
        print(f"  Evidence: {m.evidence_detail}")
        print()

    # Show sample broadMatch mappings (where no parent was found)
    print("=== SAMPLE broadMatch MAPPINGS (no suitable parent found) ===")
    broad_matches = [m for m in high_conf if m.predicate_id == "skos:broadMatch"][:3]
    for m in broad_matches:
        print(f"{m.rhea_id} -> {m.go_id}")
        print(f"  RHEA: {m.rhea_label[:70]}...")
        print(f"  GO: {m.go_label}")
        print(f"  Predicate: {m.predicate_id}")
        print(f"  Evidence: {m.evidence_detail}")
        print()

    # Show uncovered reactions
    if uncovered:
        print("=== STILL UNMAPPED ===")
        for r in uncovered:
            print(f"{r['rhea_id']}: {r['label'][:60]}...")
            print(f"  EC: {r['ec_numbers']}")
            transported = [t['name'] for t in r['transported']]
            print(f"  Transported: {transported}")
            print()


if __name__ == "__main__":
    main()
