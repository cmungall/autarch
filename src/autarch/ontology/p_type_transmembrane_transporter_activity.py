"""P-type transmembrane transporter activity.

Primary active transporter that auto-phosphorylates at a conserved Asp residue
during transport.
"""

from autarch.ontology.atpase_coupled_monoatomic_cation_transmembrane_transporter_activity import (
    ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity,
)
from autarch.ontology.go_aggregate import ExplicitGoAggregate


class PTypeTransmembraneTransporterActivity(ExplicitGoAggregate):
    """P-type transmembrane transporter activity.

    Primary active transporter that auto-phosphorylates at a conserved Asp
    residue during transport.
    """

    GO_ID = "GO:0140358"
    CONCEPT_PHRASE = "P-type transmembrane transporter activity"
    CHILD_CLASSES = (ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity,)
