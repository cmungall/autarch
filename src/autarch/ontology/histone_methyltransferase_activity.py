"""histone methyltransferase activity.

Catalysis of the reaction: S-adenosyl-L-methionine + histone = S-adenosyl-L-homocysteine + methyl-histone. Histone methylation generally occurs on either an arginine or a lysine residue.
"""

from autarch.ontology.protein_lysine_n_methyltransferase_activity import (
    _ProteinLysineMethyltransferaseBase,
)


class HistoneMethyltransferaseActivity(_ProteinLysineMethyltransferaseBase):
    """histone methyltransferase activity.

    Catalysis of the reaction: S-adenosyl-L-methionine + histone = S-adenosyl-L-homocysteine + methyl-histone. Histone methylation generally occurs on either an arginine or a lysine residue.

    The current structured RHEA cache materializes the benchmark positives as a
    generic protein-lysine scaffold rather than a histone-specific polymer.
    This class therefore detects the supported histone-methylation chemistry
    through that benchmark-aligned lysyl-protein representation.
    """

    GO_ID = "GO:0042054"
    CONCEPT_PHRASE = "histone methyltransferase activity"
