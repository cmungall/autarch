"""hydrolase acting on halide bonds in C-halide compounds.

Catalysis of the hydrolysis of any acid halide bond in substances containing
halogen atoms in organic linkage.
"""

from autarch.ontology.hydrolase_acting_on_acid_halide_bonds_in_c_halide_compounds import (
    HydrolaseActingOnAcidHalideBondsInCHalideCompounds,
)


class HydrolaseActingOnHalideBondsInCHalideCompounds(
    HydrolaseActingOnAcidHalideBondsInCHalideCompounds
):
    """hydrolase acting on halide bonds in C-halide compounds.

    Catalysis of the hydrolysis of any acid halide bond in substances
    containing halogen atoms in organic linkage.
    """

    GO_ID = "GO:0019120"
    EC_NUMBER_PREFIX = "3.8.1.-"

    def _implementation_only(self) -> None:
        """Concrete benchmark wrapper."""
