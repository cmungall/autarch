"""hydrolase acting on acid anhydrides.

Catalysis of the hydrolysis of any acid anhydride.
"""

from autarch.datamodel import Reaction
from autarch.ontology.atp_diphosphatase import ATPDiphosphatase
from autarch.ontology.atp_hydrolysis import ATPHydrolysis
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.helicase import Helicase
from autarch.ontology.hydrolase_acting_on_acid_anhydrides_in_phosphorus_containing_anhydrides import (
    HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides,
)
from autarch.ontology.nucleotide_diphosphatase import NucleotideDiphosphatase
from autarch.ontology.pyrophosphatase import Pyrophosphatase
from autarch.ontology.reaction import ReactionClass


class HydrolaseActingOnAcidAnhydrides(ReactionClass):
    """hydrolase acting on acid anhydrides.

    Catalysis of the hydrolysis of any acid anhydride.
    """

    GO_ID = "GO:0016817"
    EC_NUMBER_PREFIX = "3.6.-.-"
    CHILD_CLASSES = (
        HydrolaseActingOnAcidAnhydridesInPhosphorusContainingAnhydrides,
        NucleotideDiphosphatase,
        ATPDiphosphatase,
        Pyrophosphatase,
        ATPHydrolysis,
        Helicase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Hydrolase acting on acid anhydrides",
        )
