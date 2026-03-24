"""Mutase reaction classification using pattern DSL.

Mutases are isomerases that catalyze intramolecular group transfers.
EC 5.4.x.x classification.
"""

import re

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_COA,
    CHEBI_ATP,
    CHEBI_FAD,
    CHEBI_FADH2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    h_plus,
)


class IntramolecularPhosphotransferase(Isomerase):
    """intramolecular phosphotransferase

    Examples:
    - Phosphoglucomutase: glucose-1-P ⇌ glucose-6-P
    - Phosphoglycerate mutase: 3-phosphoglycerate ⇌ 2-phosphoglycerate
    - Methylmalonyl-CoA mutase: methylmalonyl-CoA → succinyl-CoA
    - Chorismate mutase: chorismate → prephenate
    """

    GO_ID = "GO:0016868"  # intramolecular transferase activity
    EC_NUMBER_PREFIX = "5.4.-.-"  # Intramolecular transferases (mutases)

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Simple rearrangement: substrate → rearranged substrate
        var("substrate") >> var("rearranged_product") + optional(h_plus),
        # Phosphate group migration: phospho-substrate → different-phospho-substrate
        var("phospho_substrate") >> var("rearranged_phospho_product"),
        # CoA derivative rearrangement: CoA-substrate → rearranged-CoA-product
        var("coa_substrate") >> var("rearranged_coa_product") + optional(h_plus),
    ]
    PHOSPHATE_POSITION_RE = re.compile(r"(\d+)-phosph(?:ate|o)")
    TOKEN_RE = re.compile(r"[a-z]{4,}")
    STOP_TOKENS = {
        "alpha",
        "beta",
        "phosphate",
        "phospho",
        "phosphono",
        "compound",
        "acid",
        "cis",
        "trans",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a mutase using pattern matching.

        Strategy:
        1. Exclude obvious non-isomerase chemistry (redox/ATP/CoA)
        2. Ignore balancing protons/water for stoichiometry checks
        3. Require phosphorylated substrates/products
        4. Prefer clear phosphate-position migration patterns
        """
        participants = reaction.left_participants + reaction.right_participants

        has_redox = any(
            p.chebi_id in {
                CHEBI_NAD_PLUS,
                CHEBI_NADH,
                CHEBI_NADP_PLUS,
                CHEBI_NADPH,
                CHEBI_FAD,
                CHEBI_FADH2,
                CHEBI_O2,
            }
            for p in participants
        )
        if has_redox:
            return ClassificationResult(
                is_member=False,
                explanation="Redox cofactors present - not intramolecular phosphotransferase",
            )

        has_energy_or_carrier = any(
            p.chebi_id in {CHEBI_ATP, CHEBI_COA}
            for p in participants
        )
        if has_energy_or_carrier:
            return ClassificationResult(
                is_member=False,
                explanation="ATP/CoA involvement - likely transferase, not mutase",
            )

        has_water = any(p.chebi_id == CHEBI_H2O for p in participants)
        if has_water:
            return ClassificationResult(
                is_member=False,
                explanation="Water involvement - hydration/hydrolysis rather than phosphotransfer mutase",
            )

        left_core = [
            p for p in reaction.left_participants
            if p.chebi_id not in {CHEBI_H_PLUS, CHEBI_H2O}
        ]
        right_core = [
            p for p in reaction.right_participants
            if p.chebi_id not in {CHEBI_H_PLUS, CHEBI_H2O}
        ]

        if not left_core or not right_core:
            return ClassificationResult(
                is_member=False,
                explanation="No core substrates/products after removing H+/H2O",
            )

        if len(left_core) != len(right_core) or len(left_core) > 2:
            return ClassificationResult(
                is_member=False,
                explanation="Not an intramolecular phosphotransfer stoichiometry pattern",
            )

        label = reaction.label.lower() if reaction.label else ""
        label_parts = label.split("=") if label else []
        left_label = label_parts[0] if label_parts else ""
        right_label = label_parts[1] if len(label_parts) > 1 else ""

        left_has_phosphoryl = any(p.is_phosphorylated() for p in left_core) or "phosph" in left_label
        right_has_phosphoryl = any(p.is_phosphorylated() for p in right_core) or "phosph" in right_label

        if not (left_has_phosphoryl and right_has_phosphoryl):
            return ClassificationResult(
                is_member=False,
                explanation="No phosphoryl group on both sides - not phosphotransfer mutase",
            )

        left_text = left_label
        right_text = right_label
        if not self._has_shared_core_token(left_text, right_text):
            return ClassificationResult(
                is_member=False,
                explanation="Substrate/product cores differ (no shared scaffold token)",
            )

        left_positions = set(self.PHOSPHATE_POSITION_RE.findall(left_label))
        right_positions = set(self.PHOSPHATE_POSITION_RE.findall(right_label))
        has_position_shift = bool(left_positions and right_positions and left_positions != right_positions)
        has_special_mutase_signature = (
            ("phosphoenolpyruvate" in label and "phosphonopyruvate" in label)
            or ("glyceroyl phosphate" in label and "bisphosphoglycerate" in label)
        )

        is_stereochemical_split = (
            len(left_core) == 2
            and len(right_core) == 2
            and len({p.chebi_id for p in left_core}) == 1
            and len({p.chebi_id for p in right_core}) >= 1
        )

        if len(left_core) == 1 and not (has_position_shift or has_special_mutase_signature):
            return ClassificationResult(
                is_member=False,
                explanation="No phosphate-position transfer signature for single-substrate mutase",
            )

        if len(left_core) == 2 and not (is_stereochemical_split and ("heptose" in label or has_position_shift)):
            return ClassificationResult(
                is_member=False,
                explanation="No intramolecular phosphotransfer pattern detected",
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Mutase: intramolecular phosphotransfer"
        if has_position_shift:
            explanation += " (phosphate-position migration)"
        elif is_stereochemical_split:
            explanation += " (phosphorylated stereochemical split)"
        elif has_special_mutase_signature:
            explanation += " (known phosphomutase signature)"
        else:
            explanation += " (phosphorylated intramolecular rearrangement)"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)

    @classmethod
    def _normalize_token(cls, token: str) -> str:
        for prefix in (
            "deoxy",
            "keto",
            "aldehydo",
            "mono",
            "bis",
            "phosphonato",
            "phosphono",
            "phospho",
        ):
            if token.startswith(prefix) and len(token) - len(prefix) >= 4:
                token = token[len(prefix):]
        for suffix, repl in (
            ("opyranose", "ose"),
            ("ofuranose", "ose"),
            ("heptulose", "heptose"),
            ("glyceroyl", "glycer"),
            ("glycerate", "glycer"),
            ("enolpyruvate", "pyruvate"),
        ):
            if token.endswith(suffix):
                token = token[: -len(suffix)] + repl
        return token

    @classmethod
    def _extract_core_tokens(cls, text: str) -> set[str]:
        tokens: set[str] = set()
        for raw in cls.TOKEN_RE.findall(text):
            token = cls._normalize_token(raw)
            if len(token) < 4 or token in cls.STOP_TOKENS:
                continue
            tokens.add(token)
        return tokens

    @classmethod
    def _has_shared_core_token(cls, left_text: str, right_text: str) -> bool:
        left_tokens = cls._extract_core_tokens(left_text)
        right_tokens = cls._extract_core_tokens(right_text)
        if left_tokens & right_tokens:
            return True
        for left in left_tokens:
            for right in right_tokens:
                if len(left) >= 5 and left in right:
                    return True
                if len(right) >= 5 and right in left:
                    return True
        return False
