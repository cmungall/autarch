"""N-methyltransferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.caffeine_synthase import CaffeineSynthase
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.l_histidine_n_alpha_methyltransferase import LHistidineNAlphaMethyltransferase
from autarch.ontology.nicotinamide_n_methyltransferase import NicotinamideNMethyltransferase
from autarch.ontology.phosphoethanolamine_n_methyltransferase import PhosphoethanolamineNMethyltransferase
from autarch.ontology.reaction import ReactionClass
from autarch.ontology.rs1_benzyl_1234_tetrahydroisoquinoline_n_methyltransferase import (
    RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase,
)
from autarch.ontology.trimethylsulfonium_tetrahydrofolate_n_methyltransferase import (
    TrimethylsulfoniumTetrahydrofolateNMethyltransferase,
)
from autarch.ontology.tyramine_n_methyltransferase import TyramineNMethyltransferase


class NMethyltransferase(ReactionClass):
    """N-methyltransferase activity."""

    GO_ID = "GO:0008170"
    CHILD_CLASSES = (
        NicotinamideNMethyltransferase,
        PhosphoethanolamineNMethyltransferase,
        TrimethylsulfoniumTetrahydrofolateNMethyltransferase,
        TyramineNMethyltransferase,
        RS1Benzyl1234TetrahydroisoquinolineNMethyltransferase,
        LHistidineNAlphaMethyltransferase,
        CaffeineSynthase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "N-methyltransferase",
        )
