"""N(1)-acetylpolyamine oxidase.

Catalysis of the reaction: H2O + N(1)-acetylspermine + O2 =
3-acetamidopropanal + H2O2 + spermidine.
"""

from autarch.ontology.n1_acetylpolyamine_oxidase_3_acetamidopropanal_forming import (
    N1AcetylpolyamineOxidase3AcetamidopropanalForming,
)


class N1AcetylpolyamineOxidase(N1AcetylpolyamineOxidase3AcetamidopropanalForming):
    """N(1)-acetylpolyamine oxidase.

    Catalysis of the reaction: H2O + N(1)-acetylspermine + O2 =
    3-acetamidopropanal + H2O2 + spermidine.
    """

    GO_ID = "GO:0052903"
    EC_NUMBER_PREFIX = "1.5.3.13"

    def _implementation_only(self) -> None:
        """Concrete benchmark wrapper."""
