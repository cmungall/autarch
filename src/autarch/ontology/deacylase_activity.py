"""deacylase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.d_aminoacyl_trna_deacylase import DAminoacylTRNADeacylase


class DeacylaseActivity(ExplicitGoAggregate):
    """deacylase activity."""

    GO_ID = "GO:0160215"
    EC_BROAD_XREFS = ["3.5.1.-", "3.1.1.-", "2.3.1.-"]
    CONCEPT_PHRASE = "deacylase activity"
    CHILD_CLASSES = (DAminoacylTRNADeacylase,)
