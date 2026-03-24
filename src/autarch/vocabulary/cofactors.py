"""Enzyme cofactors and coenzymes.

This module contains molecules that serve as cofactors for enzymatic reactions,
including electron carriers, group transfer agents, and prosthetic groups.
"""

from autarch.formula import molecule

# NAD/NADP system (electron carriers)
NAD = molecule("CHEBI:57540", name="NAD+")
NADH = molecule("CHEBI:57945", name="NADH")
NADP = molecule("CHEBI:58349", name="NADP+")
NADPH = molecule("CHEBI:57783", name="NADPH")

# FAD system (electron carriers)
FAD = molecule("CHEBI:57692", name="FAD")
FADH2 = molecule("CHEBI:57618", name="FADH2")
FMN = molecule("CHEBI:58210", name="FMN")
FMNH2 = molecule("CHEBI:57986", name="FMNH2")

# Coenzyme A and derivatives
CoA = molecule("CHEBI:15346", name="CoA")
acetyl_CoA = molecule("CHEBI:15351", name="acetyl-CoA")
succinyl_CoA = molecule("CHEBI:15380", name="succinyl-CoA")
malonyl_CoA = molecule("CHEBI:15531", name="malonyl-CoA")
propionyl_CoA = molecule("CHEBI:15539", name="propionyl-CoA")

# Folate system (one-carbon metabolism)
THF = molecule("CHEBI:15635", name="THF")  # Tetrahydrofolate
methyl_THF = molecule("CHEBI:15641", name="5-methyl-THF")
formyl_THF = molecule("CHEBI:15637", name="10-formyl-THF")

# S-Adenosyl methionine system (methylation)
SAM = molecule("CHEBI:15414", name="SAM")  # S-adenosylmethionine
SAH = molecule("CHEBI:16680", name="SAH")  # S-adenosylhomocysteine

# Biotin (carboxylation)
biotin = molecule("CHEBI:15956", name="biotin")
biocytin = molecule("CHEBI:16238", name="biocytin")

# Thiamine (decarboxylation)
TPP = molecule("CHEBI:9532", name="TPP")  # Thiamine pyrophosphate

# Pyridoxal phosphate (amino acid metabolism)
PLP = molecule("CHEBI:18405", name="PLP")  # Pyridoxal phosphate
PMP = molecule("CHEBI:18335", name="PMP")  # Pyridoxamine phosphate

# Lipoic acid (oxidative decarboxylation)
lipoic_acid = molecule("CHEBI:30314", name="lipoic-acid")
dihydrolipoic_acid = molecule("CHEBI:18047", name="dihydrolipoic-acid")

# Heme groups
heme = molecule("CHEBI:30413", name="heme")
heme_a = molecule("CHEBI:24479", name="heme-a")
heme_b = molecule("CHEBI:60344", name="heme-b")
heme_c = molecule("CHEBI:60562", name="heme-c")

# Metal cofactors (as complexes)
molybdopterin = molecule("CHEBI:25372", name="molybdopterin")

# Ubiquinone/Coenzyme Q
ubiquinone = molecule("CHEBI:16389", name="ubiquinone")
ubiquinol = molecule("CHEBI:17976", name="ubiquinol")

__all__ = [
    "NAD",
    "NADH",
    "NADP",
    "NADPH",
    "FAD",
    "FADH2",
    "FMN",
    "FMNH2",
    "CoA",
    "acetyl_CoA",
    "succinyl_CoA",
    "malonyl_CoA",
    "propionyl_CoA",
    "THF",
    "methyl_THF",
    "formyl_THF",
    "SAM",
    "SAH",
    "biotin",
    "biocytin",
    "TPP",
    "PLP",
    "PMP",
    "lipoic_acid",
    "dihydrolipoic_acid",
    "heme",
    "heme_a",
    "heme_b",
    "heme_c",
    "molybdopterin",
    "ubiquinone",
    "ubiquinol",
]
