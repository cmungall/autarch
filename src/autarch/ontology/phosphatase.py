"""Phosphatase reaction classification using pattern DSL.

Phosphatases are hydrolases that specifically remove phosphate groups.

History
-------

## 2025-12-22

Fixed phosphate ChEBI ID detection: RHEA uses CHEBI:16838 ("polyphosphate")
for orthophosphate products, not CHEBI:43474. Added correct ChEBI IDs.

## 2025-12-21

Refactored to use SMARTS-based phosphate detection (Moiety.PHOSPHATE) via
Participant.is_phosphorylated() instead of name-based matching for substrate
detection.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.hydrolase import Hydrolase
from autarch.ontology.reaction_diff import ReactionDiff
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_PHOSPHATE,
    adp,
    amp,
    atp,
    diphosphate,
    h_plus,
    p,
    phosphate,
    water,
)


class Phosphatase(Hydrolase):
    """phosphatase

    Examples:
    - Protein phosphatases: phospho-protein + H2O → protein + phosphate
    - Glucose-6-phosphatase: glucose-6-P + H2O → glucose + Pi
    - ATPases: ATP + H2O → ADP + Pi (energy-transducing)
    - Phospholipases: when acting on phosphate groups
    """

    GO_ID = "GO:0016791"  # phosphatase activity
    EC_NUMBER_PREFIX = "3.1.3.-"

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Generic phosphate removal: phospho-substrate + H2O → substrate + Pi
        var("phospho_substrate") + p(water)
        >> var("substrate") + p(phosphate) + optional(h_plus),
        # Diphosphate removal: substrate + H2O → product + PPi
        var("substrate") + p(water)
        >> var("product") + p(diphosphate) + optional(h_plus),
        # ATP hydrolysis (ATPase, also a phosphatase): ATP + H2O → ADP + Pi
        p(atp) + p(water) >> p(adp) + p(phosphate) + optional(h_plus),
        # ADP hydrolysis: ADP + H2O → AMP + Pi
        p(adp) + p(water) >> p(amp) + p(phosphate) + optional(h_plus),
        # ATP to AMP: ATP + H2O → AMP + PPi
        p(atp) + p(water) >> p(amp) + p(diphosphate) + optional(h_plus),
        # Multiple substrates with phosphate removal
        var("substrate1") + var("substrate2") + p(water)
        >> var("product1") + var("product2") + p(phosphate) + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a phosphatase using pattern matching.

        Strategy:
        1. First check if it's a hydrolase (parent class)
        2. Must have water as reactant (hydrolysis requirement)
        3. Must produce free phosphate/diphosphate (not just transfer)
        4. Try pattern matching for specific phosphatase patterns
        5. Verify phosphorus-containing substrate on left side
        """
        # First check if it's a hydrolase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a hydrolase: {parent_result.explanation}",
            )

        # Must have water as reactant
        water_in_reactants = any(
            p.chebi_id == CHEBI_H2O for p in reaction.left_participants
        )
        if not water_in_reactants:
            return ClassificationResult(
                is_member=False,
                explanation="No water consumption - not hydrolytic phosphate removal",
            )

        # Check for free phosphate/diphosphate production (RIGHT side only)
        # Note: RHEA uses CHEBI:16838 ("polyphosphate") for orthophosphate products
        phosphate_chebis = {
            CHEBI_PHOSPHATE,  # phosphate (CHEBI:43474)
            "CHEBI:16838",  # polyphosphate - RHEA uses this for orthophosphate!
            "CHEBI:18367",  # phosphate(3-)
        }

        diphosphate_chebis = {
            CHEBI_DIPHOSPHATE,  # diphosphate (CHEBI:33019)
            "CHEBI:18036",  # triphosphate(5-)
        }

        produces_phosphate = any(
            p.chebi_id in phosphate_chebis for p in reaction.right_participants
        )

        produces_diphosphate = any(
            p.chebi_id in diphosphate_chebis for p in reaction.right_participants
        )

        if not produces_phosphate and not produces_diphosphate:
            return ClassificationResult(
                is_member=False,
                explanation="No free phosphate/diphosphate produced - not dephosphorylation",
            )

        # CRITICAL EXCLUSION: Phosphoric anhydride hydrolases (EC 3.6.x.x) are NOT phosphatases
        # Phosphatases (EC 3.1.3.x) cleave C-O-P bonds in phosphomonoesters
        # NTPases/ATPases cleave P-O-P bonds in phosphoric anhydrides
        # Pattern: NTP + H2O → NDP + Pi (cleaves terminal P-O-P bond)
        # Pattern: NTP + H2O → NMP + PPi (cleaves both terminal P-O-P bonds)

        # Nucleoside triphosphates (substrates for P-O-P hydrolysis)
        NTP_CHEBIS = {
            CHEBI_ATP,  # ATP
            CHEBI_GTP,  # GTP
            "CHEBI:15422",  # ATP (alternate)
            "CHEBI:37565",  # UTP
            "CHEBI:17677",  # CTP
            "CHEBI:61481",  # dATP
            "CHEBI:61429",  # dGTP
            "CHEBI:61555",  # dCTP
            "CHEBI:61560",  # dTTP
            "CHEBI:60471",  # dUTP
        }

        # Nucleoside diphosphates (products of NTP hydrolysis)
        NDP_CHEBIS = {
            CHEBI_ADP,  # ADP
            CHEBI_GDP,  # GDP
            "CHEBI:17659",  # UDP
            "CHEBI:17239",  # CDP
            "CHEBI:58245",  # dADP
            "CHEBI:58305",  # dGDP
            "CHEBI:58593",  # dCDP
            "CHEBI:58369",  # dTDP
            "CHEBI:60471",  # dUDP
        }

        # Nucleoside monophosphates (products of NTP → NMP + PPi)
        NMP_CHEBIS = {
            CHEBI_AMP,  # AMP
            "CHEBI:17345",  # GMP
            "CHEBI:16695",  # UMP
            "CHEBI:17361",  # CMP
            "CHEBI:58245",  # dAMP
            "CHEBI:58115",  # dGMP
            "CHEBI:85643",  # dCMP
            "CHEBI:63528",  # dTMP
            "CHEBI:246422",  # dUMP
        }

        has_ntp = any(p.chebi_id in NTP_CHEBIS for p in reaction.left_participants)
        has_ndp = any(p.chebi_id in NDP_CHEBIS for p in reaction.right_participants)
        has_nmp = any(p.chebi_id in NMP_CHEBIS for p in reaction.right_participants)

        # NTP + H2O → NDP + Pi is phosphoric anhydride hydrolysis, not phosphatase
        if has_ntp and has_ndp and produces_phosphate:
            return ClassificationResult(
                is_member=False,
                explanation="NTP → NDP + Pi (phosphoric anhydride hydrolase EC 3.6.x.x, not phosphatase EC 3.1.3.x)",
            )

        # NTP + H2O → NMP + PPi is also phosphoric anhydride hydrolysis
        if has_ntp and has_nmp and produces_diphosphate:
            return ClassificationResult(
                is_member=False,
                explanation="NTP → NMP + PPi (phosphoric anhydride hydrolase EC 3.6.x.x, not phosphatase EC 3.1.3.x)",
            )

        # Must NOT have phosphate in reactants (it should be produced, not consumed)
        consumes_phosphate = any(
            p.chebi_id in phosphate_chebis for p in reaction.left_participants
        )

        if consumes_phosphate:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphate consumed rather than produced - not dephosphorylation",
            )

        # Get reaction diff for additional analysis
        diff = ReactionDiff(reaction)

        # Check for phosphorus-containing molecules in reactants
        # This includes ATP, ADP, AMP, and other phosphorylated compounds
        phosphorus_containing_chebis = {
            CHEBI_ATP,  # ATP
            CHEBI_ADP,  # ADP
            CHEBI_AMP,  # AMP
            CHEBI_GTP,  # GTP
            CHEBI_GDP,  # GDP
            "CHEBI:4170",  # glucose-6-phosphate
            # Add more as needed
        }

        has_phosphorus_substrate = False

        # Check if any reactant contains phosphorus via ChEBI ID or SMARTS
        for participant in reaction.left_participants:
            if participant.chebi_id in phosphorus_containing_chebis:
                has_phosphorus_substrate = True
                break
            # Use SMARTS-based phosphate detection (structural, not name-based)
            if participant.is_phosphorylated():
                has_phosphorus_substrate = True
                break

        # Also try element detection if molecules have formulas
        if not has_phosphorus_substrate and "P" in diff.reactant_elements:
            has_phosphorus_substrate = True

        if not has_phosphorus_substrate:
            return ClassificationResult(
                is_member=False, explanation="No phosphorylated substrate in reactants"
            )

        # Apply exclusions for non-phosphatase reactions that also produce phosphate
        
        # Exclude transport reactions (ATP-powered ion transport)
        # Skip ion transport detection - would need specific ChEBI IDs for ions
        
        # Exclude ligase reactions (ATP-powered bond formation)
        # Pattern: multiple substrates + ATP → product + ADP + Pi
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        has_adp = any(p.chebi_id == CHEBI_ADP for p in reaction.right_participants)
        
        if has_atp and has_adp:
            # Count non-ATP, non-water substrates
            non_atp_water_reactants = [
                p for p in reaction.left_participants
                if p.chebi_id not in {CHEBI_ATP, CHEBI_H2O}
            ]
            
            # If there are 2+ other substrates, likely a ligase reaction
            if len(non_atp_water_reactants) >= 2:
                return ClassificationResult(
                    is_member=False,
                    explanation="ATP-powered synthesis - ligase not phosphatase"
                )
            
            # Also check for patterns suggesting synthesis rather than hydrolysis
            if len(non_atp_water_reactants) >= 1:
                substrate = non_atp_water_reactants[0]
                # Look for synthesis patterns
                if substrate.name and any(pattern in substrate.name.lower() for pattern in [
                    "proline", "pyruvate", "oxaloacetate", "methylhydantoin"  # common synthesis substrates
                ]):
                    return ClassificationResult(
                        is_member=False,
                        explanation="ATP-assisted synthesis - not simple phosphate hydrolysis"
                    )
                
                # Synthesis product detection would need specific ChEBI IDs

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        # Build explanation based on what we found
        extra_info = []

        # Check if it's ATP/ADP hydrolysis (energy-transducing)
        atp_chebis = {CHEBI_ATP, CHEBI_ADP, CHEBI_AMP}
        has_atp = any(p.chebi_id in atp_chebis for p in reaction.left_participants)

        if has_atp:
            if CHEBI_ATP in [p.chebi_id for p in reaction.left_participants]:
                if CHEBI_ADP in [p.chebi_id for p in reaction.right_participants]:
                    extra_info.append("ATP→ADP")
                elif CHEBI_AMP in [
                    p.chebi_id for p in reaction.right_participants
                ]:
                    extra_info.append("ATP→AMP")
            elif CHEBI_ADP in [p.chebi_id for p in reaction.left_participants]:
                extra_info.append("ADP→AMP")

        # Check if it's protein dephosphorylation (has nitrogen)
        elif "N" in diff.reactant_elements:
            extra_info.append("likely protein substrate")

        # Check fragmentation
        if diff.is_fragmentation:
            extra_info.append("fragmentation")

        explanation = f"Phosphatase: phosphate removal ({diff.n_reactant_molecules}→{diff.n_product_molecules})"
        if extra_info:
            explanation += f", {'/'.join(extra_info)}"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
