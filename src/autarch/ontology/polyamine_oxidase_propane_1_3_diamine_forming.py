"""polyamine oxidase (propane-1,3-diamine-forming).

Catalysis of the reaction: H2O + O2 + polyamine =
aldehyde + propane-1,3-diamine + H2O2.
"""

from autarch.ontology.spermine_oxidase_propane_1_3_diamine_forming import (
    SpermineOxidasePropane13DiamineForming,
)


class PolyamineOxidasePropane13DiamineForming(
    SpermineOxidasePropane13DiamineForming
):
    """polyamine oxidase (propane-1,3-diamine-forming).

    Catalysis of the reaction: H2O + O2 + polyamine =
    aldehyde + propane-1,3-diamine + H2O2.
    """

    GO_ID = "GO:0052900"
    EC_NUMBER_PREFIX = "1.5.3.14"

    def _implementation_only(self) -> None:
        """Concrete benchmark wrapper."""
