"""S-methyltransferase activity."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.five_methyltetrahydropteroyltriglutamate_homocysteine_s_methyltransferase import FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase
from autarch.ontology.methionine_s_methyltransferase import MethionineSMethyltransferase


class SMethyltransferase(ExplicitGoAggregate):
    """S-methyltransferase activity."""

    GO_ID = "GO:0008172"
    CONCEPT_PHRASE = "S-methyltransferase activity"
    CHILD_CLASSES = (MethionineSMethyltransferase, FiveMethyltetrahydropteroyltriglutamateHomocysteineSMethyltransferase,)
