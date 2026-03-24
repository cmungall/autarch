"""Central metabolites and biomolecules.

This module contains key metabolites from central carbon metabolism,
amino acids, and other important biomolecules.
"""

from autarch.formula import molecule

# Sugars
glucose = molecule("CHEBI:17234", name="glucose")
fructose = molecule("CHEBI:28757", name="fructose")
galactose = molecule("CHEBI:28260", name="galactose")
ribose = molecule("CHEBI:33942", name="ribose")
xylose = molecule("CHEBI:17057", name="xylose")
mannose = molecule("CHEBI:4208", name="mannose")

# Sugar phosphates
glucose_6P = molecule("CHEBI:14314", name="G6P")
fructose_6P = molecule("CHEBI:15946", name="F6P")
fructose_16BP = molecule("CHEBI:16905", name="F-1,6-BP")
glucose_1P = molecule("CHEBI:16077", name="G1P")
ribose_5P = molecule("CHEBI:17797", name="R5P")

# Glycolysis intermediates
glyceraldehyde_3P = molecule("CHEBI:17138", name="G3P")
dihydroxyacetone_P = molecule("CHEBI:16108", name="DHAP")
glycerate_13BP = molecule("CHEBI:16001", name="1,3-BPG")
glycerate_3P = molecule("CHEBI:17794", name="3PG")
glycerate_2P = molecule("CHEBI:17835", name="2PG")
phosphoenolpyruvate = molecule("CHEBI:18021", name="PEP")

# Organic acids
pyruvate = molecule("CHEBI:15361", name="pyruvate")
lactate = molecule("CHEBI:24996", name="lactate")
acetate = molecule("CHEBI:30089", name="acetate")
formate = molecule("CHEBI:15740", name="formate")

# TCA cycle intermediates
citrate = molecule("CHEBI:16947", name="citrate")
isocitrate = molecule("CHEBI:16219", name="isocitrate")
alpha_ketoglutarate = molecule("CHEBI:16810", name="α-KG")
succinate = molecule("CHEBI:15741", name="succinate")
fumarate = molecule("CHEBI:18012", name="fumarate")
malate = molecule("CHEBI:15589", name="malate")
oxaloacetate = molecule("CHEBI:16452", name="OAA")

# Amino acids (20 standard)
alanine = molecule("CHEBI:16449", name="alanine")
arginine = molecule("CHEBI:16467", name="arginine")
asparagine = molecule("CHEBI:17196", name="asparagine")
aspartate = molecule("CHEBI:35391", name="aspartate")
cysteine = molecule("CHEBI:17561", name="cysteine")
glutamate = molecule("CHEBI:16015", name="glutamate")
glutamine = molecule("CHEBI:18050", name="glutamine")
glycine = molecule("CHEBI:15428", name="glycine")
histidine = molecule("CHEBI:15971", name="histidine")
isoleucine = molecule("CHEBI:17191", name="isoleucine")
leucine = molecule("CHEBI:15603", name="leucine")
lysine = molecule("CHEBI:18019", name="lysine")
methionine = molecule("CHEBI:16643", name="methionine")
phenylalanine = molecule("CHEBI:17295", name="phenylalanine")
proline = molecule("CHEBI:17203", name="proline")
serine = molecule("CHEBI:17115", name="serine")
threonine = molecule("CHEBI:16857", name="threonine")
tryptophan = molecule("CHEBI:16828", name="tryptophan")
tyrosine = molecule("CHEBI:17895", name="tyrosine")
valine = molecule("CHEBI:16414", name="valine")

# Lipids and fatty acids
palmitate = molecule("CHEBI:15756", name="palmitate")
stearate = molecule("CHEBI:28842", name="stearate")
oleate = molecule("CHEBI:16196", name="oleate")
glycerol = molecule("CHEBI:17754", name="glycerol")
glycerol_3P = molecule("CHEBI:15978", name="glycerol-3P")

# Nucleotides/nucleosides
adenosine = molecule("CHEBI:16335", name="adenosine")
guanosine = molecule("CHEBI:16750", name="guanosine")
cytidine = molecule("CHEBI:17562", name="cytidine")
uridine = molecule("CHEBI:16704", name="uridine")
thymidine = molecule("CHEBI:17748", name="thymidine")

# Other important metabolites
urea = molecule("CHEBI:16199", name="urea")
creatine = molecule("CHEBI:16919", name="creatine")
creatinine = molecule("CHEBI:16737", name="creatinine")

__all__ = [
    # Sugars
    "glucose",
    "fructose",
    "galactose",
    "ribose",
    "xylose",
    "mannose",
    # Sugar phosphates
    "glucose_6P",
    "fructose_6P",
    "fructose_16BP",
    "glucose_1P",
    "ribose_5P",
    # Glycolysis
    "glyceraldehyde_3P",
    "dihydroxyacetone_P",
    "glycerate_13BP",
    "glycerate_3P",
    "glycerate_2P",
    "phosphoenolpyruvate",
    # Organic acids
    "pyruvate",
    "lactate",
    "acetate",
    "formate",
    # TCA cycle
    "citrate",
    "isocitrate",
    "alpha_ketoglutarate",
    "succinate",
    "fumarate",
    "malate",
    "oxaloacetate",
    # Amino acids
    "alanine",
    "arginine",
    "asparagine",
    "aspartate",
    "cysteine",
    "glutamate",
    "glutamine",
    "glycine",
    "histidine",
    "isoleucine",
    "leucine",
    "lysine",
    "methionine",
    "phenylalanine",
    "proline",
    "serine",
    "threonine",
    "tryptophan",
    "tyrosine",
    "valine",
    # Lipids
    "palmitate",
    "stearate",
    "oleate",
    "glycerol",
    "glycerol_3P",
    # Nucleotides
    "adenosine",
    "guanosine",
    "cytidine",
    "uridine",
    "thymidine",
    # Others
    "urea",
    "creatine",
    "creatinine",
]
