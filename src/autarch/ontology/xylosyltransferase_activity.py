"""xylosyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.udp_xylosyltransferase_activity import (
    UDPXylosyltransferaseActivity,
)


class XylosyltransferaseActivity(ExplicitGoAggregate):
    """xylosyltransferase activity."""

    GO_ID = "GO:0042285"
    CONCEPT_PHRASE = "xylosyltransferase activity"
    CHILD_CLASSES = (UDPXylosyltransferaseActivity,)
