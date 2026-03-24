"""cis-trans isomerase activity."""

from autarch.ontology.cis_trans_isomerases import CisTransIsomerases
from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.peptidyl_prolyl_cis_trans_isomerase import PeptidylProlylCisTransIsomerase


class CisTransIsomerase(ExplicitGoAggregate):
    """cis-trans isomerase activity."""

    GO_ID = "GO:0016859"
    EC_NUMBER_PREFIX = "5.2.-.-"
    CHILD_CLASSES = (
        CisTransIsomerases,
        PeptidylProlylCisTransIsomerase,
    )
