"""Common molecule constants for pattern matching.

This module provides:
1. ChEBI ID constants (strings) for direct comparisons
2. Pre-defined Participant objects for pattern matching
3. Semantic groups of related molecules
4. Helper functions for participant matching

Examples:
    >>> from autarch.molecules import atp, adp, p, CHEBI_ATP, has_chebi
    >>> from autarch.pattern_dsl import var
    >>>
    >>> # Create a pattern using operator overloading
    >>> pattern = p(atp) + var("S") >> p(adp) + var("P")
    >>> len(pattern.left_participants)
    2
    >>>
    >>> # Check if a participant is ATP
    >>> has_chebi(atp, CHEBI_ATP)
    True
"""

from autarch.datamodel import Participant
from typing import TYPE_CHECKING, List, Set

# =============================================================================
# ChEBI ID Constants (bare strings for direct comparison)
# =============================================================================

# Small molecules
CHEBI_H2O = "CHEBI:15377"
CHEBI_H_PLUS = "CHEBI:15378"
CHEBI_CO2 = "CHEBI:16526"
CHEBI_O2 = "CHEBI:15379"
CHEBI_H2O2 = "CHEBI:16240"
CHEBI_NH3 = "CHEBI:16134"
CHEBI_NH4 = "CHEBI:28938"
CHEBI_PHOSPHATE = "CHEBI:43474"
CHEBI_DIPHOSPHATE = "CHEBI:33019"

# Energy carriers (nucleoside phosphates)
CHEBI_ATP = "CHEBI:30616"
CHEBI_ADP = "CHEBI:456216"
CHEBI_AMP = "CHEBI:456215"
CHEBI_GTP = "CHEBI:37565"
CHEBI_GDP = "CHEBI:58189"
CHEBI_GMP = "CHEBI:58115"
CHEBI_UTP = "CHEBI:46398"
CHEBI_UDP = "CHEBI:58223"
CHEBI_UMP = "CHEBI:57865"
CHEBI_CTP = "CHEBI:37563"
CHEBI_CDP = "CHEBI:58069"
CHEBI_CMP = "CHEBI:60377"

# Redox cofactors
CHEBI_NAD_PLUS = "CHEBI:57540"
CHEBI_NADH = "CHEBI:57945"
CHEBI_NADP_PLUS = "CHEBI:58349"
CHEBI_NADPH = "CHEBI:57783"
CHEBI_FAD = "CHEBI:57692"
CHEBI_FADH2 = "CHEBI:57618"
CHEBI_FMN = "CHEBI:58210"
CHEBI_FMNH2 = "CHEBI:57618"  # Same as FADH2 for reduced form

# Methyl/acyl transfer cofactors
CHEBI_SAM = "CHEBI:59789"
CHEBI_SAH = "CHEBI:57856"
CHEBI_COA = "CHEBI:57287"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_MALONYL_COA = "CHEBI:57384"
CHEBI_PALMITOYL_COA = "CHEBI:57379"

# =============================================================================
# Semantic Groups (sets of ChEBI IDs for category matching)
# =============================================================================

# Redox cofactor pairs (NAD+/NADH, NADP+/NADPH, FAD/FADH2)
REDOX_COFACTORS: Set[str] = {
    CHEBI_NAD_PLUS, CHEBI_NADH,
    CHEBI_NADP_PLUS, CHEBI_NADPH,
    CHEBI_FAD, CHEBI_FADH2,
    CHEBI_FMN, CHEBI_FMNH2,
}

# Cofactor pairs for detecting redox reactions
REDOX_COFACTOR_PAIRS: Set[tuple] = {
    (CHEBI_NAD_PLUS, CHEBI_NADH),
    (CHEBI_NADP_PLUS, CHEBI_NADPH),
    (CHEBI_FAD, CHEBI_FADH2),
}

# Phosphate donors (nucleoside triphosphates)
PHOSPHATE_DONORS: Set[str] = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}

# Phosphate acceptor products (nucleoside diphosphates)
PHOSPHATE_PRODUCTS: Set[str] = {CHEBI_ADP, CHEBI_GDP, CHEBI_UDP, CHEBI_CDP}

# All nucleoside phosphates (for transferase detection)
NUCLEOSIDE_PHOSPHATES: Set[str] = {
    CHEBI_ATP, CHEBI_ADP, CHEBI_AMP,
    CHEBI_GTP, CHEBI_GDP, CHEBI_GMP,
    CHEBI_UTP, CHEBI_UDP, CHEBI_UMP,
    CHEBI_CTP, CHEBI_CDP, CHEBI_CMP,
}

# Glycosyl transfer carriers (NDP-sugars release NDP)
GLYCOSYL_CARRIERS: Set[str] = {CHEBI_UDP, CHEBI_GDP, CHEBI_CDP}

# Acyl transfer carriers
ACYL_CARRIERS: Set[str] = {CHEBI_COA, CHEBI_ACETYL_COA, CHEBI_MALONYL_COA, CHEBI_PALMITOYL_COA}

# Methyl transfer pair
METHYL_TRANSFER_PAIR = (CHEBI_SAM, CHEBI_SAH)

# Small molecules often present as spectators or water-mediated
SMALL_MOLECULES: Set[str] = {CHEBI_H2O, CHEBI_H_PLUS, CHEBI_CO2, CHEBI_O2, CHEBI_NH3, CHEBI_NH4}

# =============================================================================
# Helper Functions
# =============================================================================


def has_chebi(participant: Participant, chebi_id: str) -> bool:
    """Check if participant matches a specific ChEBI ID.

    Args:
        participant: Participant to check
        chebi_id: ChEBI ID to match (e.g., "CHEBI:30616")

    Returns:
        True if participant has the specified ChEBI ID

    Examples:
        >>> from autarch.molecules import atp, CHEBI_ATP, has_chebi
        >>> has_chebi(atp, CHEBI_ATP)
        True
        >>> has_chebi(atp, CHEBI_ADP)
        False
    """
    return participant.chebi_id == chebi_id


def has_any_chebi(participant: Participant, chebi_ids: Set[str]) -> bool:
    """Check if participant matches any of the ChEBI IDs in a set.

    Args:
        participant: Participant to check
        chebi_ids: Set of ChEBI IDs to match against

    Returns:
        True if participant's ChEBI ID is in the set

    Examples:
        >>> from autarch.molecules import atp, PHOSPHATE_DONORS, has_any_chebi
        >>> has_any_chebi(atp, PHOSPHATE_DONORS)
        True
    """
    return participant.chebi_id in chebi_ids


def find_by_chebi(participants: List[Participant], chebi_id: str) -> List[Participant]:
    """Find all participants with a specific ChEBI ID.

    Args:
        participants: List of participants to search
        chebi_id: ChEBI ID to find

    Returns:
        List of matching participants (may be empty)

    Examples:
        >>> from autarch.molecules import atp, adp, find_by_chebi, CHEBI_ATP
        >>> matches = find_by_chebi([atp, adp], CHEBI_ATP)
        >>> len(matches)
        1
        >>> matches[0].chebi_id
        'CHEBI:30616'
    """
    return [p for p in participants if p.chebi_id == chebi_id]


def find_any_chebi(participants: List[Participant], chebi_ids: Set[str]) -> List[Participant]:
    """Find all participants matching any of the ChEBI IDs.

    Args:
        participants: List of participants to search
        chebi_ids: Set of ChEBI IDs to match

    Returns:
        List of matching participants

    Examples:
        >>> from autarch.molecules import atp, adp, water, find_any_chebi, PHOSPHATE_DONORS
        >>> matches = find_any_chebi([atp, adp, water], PHOSPHATE_DONORS)
        >>> len(matches)
        1
    """
    return [p for p in participants if p.chebi_id in chebi_ids]


def has_redox_cofactor(participants: List[Participant]) -> bool:
    """Check if any participant is a redox cofactor (NAD+/NADH, etc.).

    Examples:
        >>> from autarch.molecules import nad_plus, water, has_redox_cofactor
        >>> has_redox_cofactor([nad_plus, water])
        True
    """
    return any(p.chebi_id in REDOX_COFACTORS for p in participants)


def has_phosphate_donor(participants: List[Participant]) -> bool:
    """Check if any participant is a phosphate donor (ATP, GTP, etc.).

    Examples:
        >>> from autarch.molecules import atp, water, has_phosphate_donor
        >>> has_phosphate_donor([atp, water])
        True
    """
    return any(p.chebi_id in PHOSPHATE_DONORS for p in participants)


# =============================================================================
# Participant Objects (for pattern matching)
# =============================================================================

if TYPE_CHECKING:
    from autarch.pattern_dsl import PatternParticipant


def p(participant: Participant) -> "PatternParticipant":
    """Convert a Participant to PatternParticipant for pattern DSL.

    This convenience function enables operator overloading in patterns.

    Args:
        participant: Regular Participant to convert

    Returns:
        PatternParticipant for use in patterns

    Examples:
        >>> pattern_atp = p(atp)
        >>> pattern_atp.chebi_id
        'CHEBI:30616'
    """
    from autarch.pattern_dsl import to_pattern

    return to_pattern(participant)


# Energy carriers
atp = ATP = Participant(chebi_id="CHEBI:30616", name="ATP")
adp = ADP = Participant(chebi_id="CHEBI:456216", name="ADP")
amp = AMP = Participant(chebi_id="CHEBI:456215", name="AMP")
gtp = GTP = Participant(chebi_id="CHEBI:37565", name="GTP")
gdp = GDP = Participant(chebi_id="CHEBI:58189", name="GDP")
gmp = GMP = Participant(chebi_id="CHEBI:58115", name="GMP")

# Cofactors
nad_plus = NAD_PLUS = Participant(chebi_id="CHEBI:57540", name="NAD+")
nadh = NADH = Participant(chebi_id="CHEBI:57945", name="NADH")
nadp_plus = NADP_PLUS = Participant(chebi_id="CHEBI:58349", name="NADP+")
nadph = NADPH = Participant(chebi_id="CHEBI:57783", name="NADPH")
fad = FAD = Participant(chebi_id="CHEBI:57692", name="FAD")
fadh2 = FADH2 = Participant(chebi_id="CHEBI:57618", name="FADH2")
coenzyme_a = CoA = Participant(chebi_id="CHEBI:57287", name="CoA")
acetyl_coa = ACETYL_COA = Participant(chebi_id="CHEBI:57288", name="acetyl-CoA")

# Methyl transfer cofactors
sam = SAM = Participant(chebi_id="CHEBI:59789", name="S-adenosyl-L-methionine")
sah = SAH = Participant(chebi_id="CHEBI:57856", name="S-adenosyl-L-homocysteine")

# Small molecules
water = H2O = Participant(chebi_id="CHEBI:15377", name="H2O", smiles="O")
h_plus = proton = H_PLUS = Participant(chebi_id="CHEBI:15378", name="H+", smiles="[H+]")
phosphate = Pi = PHOSPHATE = Participant(
    chebi_id="CHEBI:43474", name="phosphate", smiles="[O-]P([O-])([O-])=O"
)
diphosphate = PPi = DIPHOSPHATE = Participant(
    chebi_id="CHEBI:33019", name="diphosphate", smiles="[O-]P([O-])(=O)OP([O-])([O-])=O"
)
co2 = CO2 = Participant(chebi_id="CHEBI:16526", name="CO2", smiles="O=C=O")
oxygen = O2 = Participant(chebi_id="CHEBI:15379", name="O2", smiles="O=O")
hydrogen_peroxide = H2O2 = Participant(chebi_id="CHEBI:16240", name="H2O2", smiles="OO")
ammonia = NH3 = Participant(chebi_id="CHEBI:16134", name="NH3", smiles="N")
ammonium = NH4_PLUS = Participant(chebi_id="CHEBI:28938", name="NH4+", smiles="[NH4+]")

# Common metabolites
glucose = GLUCOSE = Participant(chebi_id="CHEBI:17234", name="glucose")
glucose_6_phosphate = G6P = Participant(
    chebi_id="CHEBI:17665", name="glucose-6-phosphate"
)
fructose_6_phosphate = F6P = Participant(
    chebi_id="CHEBI:57579", name="fructose-6-phosphate"
)
pyruvate = PYRUVATE = Participant(chebi_id="CHEBI:15361", name="pyruvate")
lactate = LACTATE = Participant(chebi_id="CHEBI:24996", name="lactate")
acetate = ACETATE = Participant(chebi_id="CHEBI:30089", name="acetate")
citrate = CITRATE = Participant(chebi_id="CHEBI:30769", name="citrate")
alpha_ketoglutarate = AKG = Participant(chebi_id="CHEBI:16810", name="2-oxoglutarate")
succinate = SUCCINATE = Participant(chebi_id="CHEBI:30031", name="succinate")
malate = MALATE = Participant(chebi_id="CHEBI:30797", name="malate")
oxaloacetate = OAA = Participant(chebi_id="CHEBI:30744", name="oxaloacetate")

# Amino acids (common ones)
glutamate = GLU = Participant(chebi_id="CHEBI:29985", name="glutamate")
glutamine = GLN = Participant(chebi_id="CHEBI:29987", name="glutamine")
aspartate = ASP = Participant(chebi_id="CHEBI:29991", name="aspartate")
alanine = ALA = Participant(chebi_id="CHEBI:57972", name="alanine")
glycine = GLY = Participant(chebi_id="CHEBI:57305", name="glycine")
serine = SER = Participant(chebi_id="CHEBI:33384", name="serine")

# Lipids and fatty acids
palmitate = PALMITATE = Participant(chebi_id="CHEBI:29032", name="palmitate")
oleate = OLEATE = Participant(chebi_id="CHEBI:25664", name="oleate")
cholesterol = CHOLESTEROL = Participant(chebi_id="CHEBI:16113", name="cholesterol")

# Ions
sodium_ion = NA_PLUS = Participant(chebi_id="CHEBI:29101", name="Na+")
potassium_ion = K_PLUS = Participant(chebi_id="CHEBI:29103", name="K+")
calcium_ion = CA_2PLUS = Participant(chebi_id="CHEBI:29108", name="Ca2+")
chloride_ion = CL_MINUS = Participant(chebi_id="CHEBI:17996", name="Cl-")
magnesium_ion = MG_2PLUS = Participant(chebi_id="CHEBI:18420", name="Mg2+")
iron_ion_fe2 = FE_2PLUS = Participant(chebi_id="CHEBI:29033", name="Fe2+")
iron_ion_fe3 = FE_3PLUS = Participant(chebi_id="CHEBI:29034", name="Fe3+")

# Nucleotides (besides ATP/ADP/AMP)
utp = UTP = Participant(chebi_id="CHEBI:46398", name="UTP")
udp = UDP = Participant(chebi_id="CHEBI:58223", name="UDP")
ump = UMP = Participant(chebi_id="CHEBI:57865", name="UMP")
ctp = CTP = Participant(chebi_id="CHEBI:37563", name="CTP")
cdp = CDP = Participant(chebi_id="CHEBI:58069", name="CDP")
cmp = CMP = Participant(chebi_id="CHEBI:60377", name="CMP")

# DNA/RNA building blocks
datp = dATP = Participant(chebi_id="CHEBI:61404", name="dATP")
dgtp = dGTP = Participant(chebi_id="CHEBI:61429", name="dGTP")
dctp = dCTP = Participant(chebi_id="CHEBI:61481", name="dCTP")
dttp = dTTP = Participant(chebi_id="CHEBI:63528", name="dTTP")
