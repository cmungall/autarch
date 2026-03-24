"""histone modifying activity.

A catalytic activity that acts on a histone protein. Reversible histone modifications contribute to regulation of gene expression.
"""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.histone_methyltransferase_activity import (
    HistoneMethyltransferaseActivity,
)


class HistoneModifyingActivity(ExplicitGoAggregate):
    """histone modifying activity.

    A catalytic activity that acts on a histone protein. Reversible histone modifications contribute to regulation of gene expression.
    """

    GO_ID = "GO:0140993"
    CONCEPT_PHRASE = "histone modifying activity"
    CHILD_CLASSES = (HistoneMethyltransferaseActivity,)
