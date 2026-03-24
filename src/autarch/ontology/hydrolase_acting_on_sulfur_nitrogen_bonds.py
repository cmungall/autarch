"""hydrolase acting on sulfur-nitrogen bonds.

Catalysis of the hydrolysis of acid sulfur-nitrogen bonds.
"""

from autarch.ontology.hydrolase_acting_on_acid_sulfur_nitrogen_bonds import (
    HydrolaseActingOnAcidSulfurNitrogenBonds,
)


class HydrolaseActingOnSulfurNitrogenBonds(
    HydrolaseActingOnAcidSulfurNitrogenBonds
):
    """hydrolase acting on sulfur-nitrogen bonds.

    Catalysis of the hydrolysis of acid sulfur-nitrogen bonds.
    """

    GO_ID = "GO:0016826"
    EC_NUMBER_PREFIX = "3.10.-.-"

    def _implementation_only(self) -> None:
        """Concrete benchmark wrapper."""
