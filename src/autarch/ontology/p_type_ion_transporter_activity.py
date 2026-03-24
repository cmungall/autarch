"""P-type ion transporter activity.

Enables ATP-driven ion transport via a phosphorylated Asp transporter
intermediate.
"""

from autarch.ontology.atpase_coupled_monoatomic_cation_transmembrane_transporter_activity import (
    ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity,
)
from autarch.ontology.go_aggregate import ExplicitGoAggregate


class PTypeIonTransporterActivity(ExplicitGoAggregate):
    """P-type ion transporter activity.

    Enables ATP-driven ion transport via a phosphorylated Asp transporter
    intermediate.
    """

    GO_ID = "GO:0015662"
    CONCEPT_PHRASE = "P-type ion transporter activity"
    CHILD_CLASSES = (ATPaseCoupledMonoatomicCationTransmembraneTransporterActivity,)
