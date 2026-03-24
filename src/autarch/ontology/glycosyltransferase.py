"""Glycosyltransferase reaction classification using pattern DSL.

Glycosyltransferases transfer glycosyl groups from donors to acceptors.
EC 2.4.x.x classification.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import IDENTIFIER_FRIENDLY_EVIDENCE, ReactionClass
from autarch.pattern_dsl import var, optional, match_patterns
from autarch.molecules import (
    h_plus, CHEBI_UDP, CHEBI_GDP, CHEBI_ADP, CHEBI_CDP,
    CHEBI_H2O, CHEBI_H_PLUS, CHEBI_ATP, CHEBI_AMP, CHEBI_PHOSPHATE,
    CHEBI_UTP, CHEBI_GTP, CHEBI_DIPHOSPHATE, CHEBI_CO2,
)


class Glycosyltransferase(ReactionClass):
    """glycosyltransferase

    Examples:
    - Glucosyltransferase: UDP-glucose + acceptor → UDP + glucoside
    - Galactosyltransferase: UDP-galactose + acceptor → UDP + galactoside
    - Fucosyltransferase: GDP-fucose + acceptor → GDP + fucoside
    - Glycogen synthase: UDP-glucose + glycogen → UDP + extended glycogen
    """

    GO_ID = "GO:0016757"  # transferase activity, transferring glycosyl groups
    EC_NUMBER_PREFIX = "2.4.-.-"  # Glycosyltransferases
    EVALUATION_EVIDENCE = IDENTIFIER_FRIENDLY_EVIDENCE
    PHOSPHATE_IDS = {CHEBI_PHOSPHATE, "CHEBI:16838"}
    DIPHOSPHATE_IDS = {CHEBI_DIPHOSPHATE, "CHEBI:18036"}
    CYCLIC_NUCLEOTIDES = {"CHEBI:58165", "CHEBI:57746"}  # cAMP, cGMP
    PEP_IDS = {"CHEBI:58702"}  # phosphoenolpyruvate
    NUCLEOTIDE_SUGAR_IDS = {
        # UDP-sugars
        "CHEBI:18066",
        "CHEBI:18307",
        "CHEBI:58885",
        "CHEBI:66914",
        "CHEBI:57705",
        "CHEBI:58052",
        "CHEBI:57632",
        "CHEBI:67138",
        "CHEBI:68623",
        "CHEBI:83836",
        # GDP-sugars
        "CHEBI:16264",
        "CHEBI:15819",
        "CHEBI:57527",
        "CHEBI:57273",
        "CHEBI:62230",
        # CMP/ADP/NDP sugars
        "CHEBI:57812",
        "CHEBI:57498",
        "CHEBI:61506",
        "CHEBI:76533",
        "CHEBI:57930",
    }
    FREE_NUCLEOTIDE_IDS = {
        CHEBI_UDP,
        CHEBI_GDP,
        CHEBI_ADP,
        CHEBI_UTP,
        CHEBI_CDP,
        "CHEBI:57865",  # UMP
        "CHEBI:58115",  # GMP
        "CHEBI:60377",  # CMP
        "CHEBI:58369",  # dTDP
        CHEBI_ATP,
    }

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # UDP-sugar transfer: UDP-sugar + acceptor → UDP + glycoside + H+
        var("udp_sugar") + var("acceptor") >> var("udp") + var("glycoside") + optional(h_plus),
        # GDP-sugar transfer: GDP-sugar + acceptor → GDP + glycoside + H+
        var("gdp_sugar") + var("acceptor") >> var("gdp") + var("glycoside") + optional(h_plus),
        # General nucleotide-sugar: NDP-sugar + acceptor → NDP + product + H+
        var("ndp_sugar") + var("acceptor") >> var("ndp") + var("product") + optional(h_plus),
        # Polysaccharide synthesis: NDP-sugar + polysaccharide → NDP + extended-polysaccharide
        var("sugar_donor") + var("polysaccharide") >> var("nucleotide") + var("extended_polysaccharide") + optional(h_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is a glycosyltransferase using pattern matching.

        Strategy (independent of parent Transferase):
        1. Must involve nucleotide-sugar donors (UDP-, GDP-, etc.)
        2. Must release nucleotide (UDP, GDP, etc.) as product
        3. Check for glycosyl group transfer patterns
        4. Exclude nucleotide activation/cyclase/GTPase side chemistry
        """
        # Check for nucleotide-sugar donors by ChEBI IDs.
        has_nucleotide_sugar = any(
            p.chebi_id in self.NUCLEOTIDE_SUGAR_IDS
            for p in reaction.left_participants
        )

        # Check for nucleotide products (released when sugar is transferred)
        has_nucleotide_product = any(
            p.chebi_id in self.FREE_NUCLEOTIDE_IDS
            for p in reaction.right_participants
        )

        # Must have evidence of nucleotide-sugar involvement
        if not (has_nucleotide_sugar and has_nucleotide_product):
            return ClassificationResult(
                is_member=False,
                explanation="No nucleotide-sugar transfer pattern detected"
            )

        # Apply exclusions before classifying as glycosyltransferase

        # Exclude simple hydrolysis (water + nucleotide-sugar → sugar + nucleotide)
        has_water = any(
            p.chebi_id == CHEBI_H2O
            for p in reaction.left_participants
        )

        # Count non-water, non-nucleotide substrates (must have actual acceptor molecule)
        spectator_chebi = {CHEBI_H2O, CHEBI_H_PLUS}
        non_cofactor_substrates = [
            p for p in reaction.left_participants
            if p.chebi_id not in spectator_chebi and
            p.chebi_id not in self.NUCLEOTIDE_SUGAR_IDS
        ]

        if has_water and len(non_cofactor_substrates) == 0:
            return ClassificationResult(
                is_member=False,
                explanation="Simple nucleotide-sugar hydrolysis - not glycosyltransferase"
            )

        # Exclude phosphorylase reactions (use phosphate instead of nucleotides)
        has_phosphate = any(
            p.chebi_id in self.PHOSPHATE_IDS
            for p in reaction.left_participants + reaction.right_participants
        )

        if has_phosphate and not has_nucleotide_product:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphorylase reaction - not nucleotide-dependent glycosyltransferase"
            )

        # Exclude ligase reactions (ATP consumption with bond formation)
        has_atp = any(p.chebi_id == CHEBI_ATP for p in reaction.left_participants)
        has_amp_or_adp = any(
            p.chebi_id in {CHEBI_AMP, CHEBI_ADP}for p in reaction.right_participants
        )

        if has_atp and has_amp_or_adp:
            # Could be nucleotide-sugar synthesis, not transfer
            return ClassificationResult(
                is_member=False,
                explanation="ATP-dependent ligase - likely nucleotide-sugar synthesis"
            )

        # Exclude cyclases (GTP → cGMP + PPi, ATP → cAMP + PPi)
        # These convert NTP to cyclic nucleotide, NOT glycosyl transfer
        has_diphosphate_product = any(
            p.chebi_id in self.DIPHOSPHATE_IDS
            for p in reaction.right_participants
        )
        has_cyclic_product = any(
            p.chebi_id in self.CYCLIC_NUCLEOTIDES
            for p in reaction.right_participants
        )
        if has_diphosphate_product and has_cyclic_product:
            return ClassificationResult(
                is_member=False,
                explanation="Cyclase (NTP → cyclic nucleotide + PPi) - not glycosyltransferase"
            )

        # Exclude epimerases/isomerases (UDP-sugar = UDP-epimer, no transfer)
        # These are 1:1 conversions where the sugar moiety changes configuration
        if len(reaction.left_participants) == 1 and len(reaction.right_participants) == 1:
            # Single substrate → single product is likely epimerase, not transfer
            return ClassificationResult(
                is_member=False,
                explanation="1:1 interconversion - likely epimerase not glycosyltransferase"
            )

        # Exclude PEP carboxykinase (oxaloacetate + GTP → PEP + GDP + CO2)
        # This uses GTP for phosphoryl transfer to create PEP, not glycosyl transfer
        has_co2_product = any(
            p.chebi_id in {CHEBI_CO2, "CHEBI:13283"}
            for p in reaction.right_participants
        )
        has_pep_product = any(
            p.chebi_id in self.PEP_IDS
            for p in reaction.right_participants
        )
        if has_co2_product and has_pep_product:
            return ClassificationResult(
                is_member=False,
                explanation="PEP carboxykinase (uses NTP for phosphoryl transfer) - not glycosyltransferase"
            )

        # Exclude reactions where GTP/GDP is used for energy, not sugar donation
        # Pattern: substrate + GTP → product + GDP + Pi (GTPase-coupled reaction)
        has_gtp = any(p.chebi_id == CHEBI_GTP for p in reaction.left_participants)
        has_gdp = any(p.chebi_id == CHEBI_GDP for p in reaction.right_participants)
        has_phosphate_product = any(
            p.chebi_id in self.PHOSPHATE_IDS for p in reaction.right_participants
        )
        # GTP + substrate → GDP + Pi + product is GTPase, not glycosyltransferase
        if has_gtp and has_gdp and has_phosphate_product:
            return ClassificationResult(
                is_member=False,
                explanation="GTPase-coupled reaction - not glycosyltransferase"
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Glycosyltransferase: glycosyl group transfer"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)
