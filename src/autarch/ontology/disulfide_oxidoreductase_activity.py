"""disulfide oxidoreductase activity.

Catalysis of the reaction: substrate with reduced sulfide groups = substrate
with oxidized disulfide bonds.
"""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.oxidoreductase_acting_on_a_sulfur_group_of_donors_disulfide_as_acceptor import (
    OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_ch_or_ch2_groups_disulfide_as_acceptor import (
    OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor,
)
from autarch.ontology.ribonucleoside_diphosphate_reductase_thioredoxin_disulfide_as_acceptor import (
    RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor,
)


class DisulfideOxidoreductaseActivity(ExplicitGoAggregate):
    """disulfide oxidoreductase activity.

    Catalysis of the reaction: substrate with reduced sulfide groups =
    substrate with oxidized disulfide bonds.
    """

    GO_ID = "GO:0015036"
    CONCEPT_PHRASE = "disulfide oxidoreductase activity"
    CHILD_CLASSES = (
        OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor,
        OxidoreductaseActingOnCHOrCH2GroupsDisulfideAsAcceptor,
        RibonucleosideDiphosphateReductaseThioredoxinDisulfideAsAcceptor,
    )
