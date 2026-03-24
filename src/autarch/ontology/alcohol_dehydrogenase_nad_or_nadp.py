"""alcohol dehydrogenase [NAD(P)+] activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.alcohol_dehydrogenase_nad import AlcoholDehydrogenaseNAD


class AlcoholDehydrogenaseNADOrNADP(ExplicitGoAggregate):
    """alcohol dehydrogenase [NAD(P)+] activity."""

    GO_ID = "GO:0018455"
    CONCEPT_PHRASE = "alcohol dehydrogenase [NAD(P)+] activity"
    CHILD_CLASSES = (AlcoholDehydrogenaseNAD,)
