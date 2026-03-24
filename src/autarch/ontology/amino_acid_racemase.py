"""amino-acid racemase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.racemase_and_epimerase_acting_on_amino_acids_and_derivatives import (
    RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives,
)


class AminoAcidRacemase(ExplicitGoAggregate):
    """amino-acid racemase activity."""

    GO_ID = "GO:0047661"
    CONCEPT_PHRASE = "amino-acid racemase activity"
    CHILD_CLASSES = (RacemaseAndEpimeraseActingOnAminoAcidsAndDerivatives,)
