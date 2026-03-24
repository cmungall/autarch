"""acetyltransferase activity."""

from autarch.datamodel import Reaction
from autarch.ontology.ec_prefix_aggregate import (
    aggregate_explicit_child_membership,
    aggregate_explicit_child_supports_evaluation,
)
from autarch.ontology.n_terminal_amino_acid_n_alpha_acetyltransferase_nat_a import (
    NTerminalAminoAcidNAlphaAcetyltransferaseNatA,
)
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_b import (
    NTerminalMethionineNAlphaAcetyltransferaseNatB,
)
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_c import (
    NTerminalMethionineNAlphaAcetyltransferaseNatC,
)
from autarch.ontology.n_terminal_methionine_n_alpha_acetyltransferase_nat_e import (
    NTerminalMethionineNAlphaAcetyltransferaseNatE,
)
from autarch.ontology.polysialic_acid_o_acetyltransferase import (
    PolysialicAcidOAcetyltransferase,
)
from autarch.ontology.reaction import ReactionClass


class Acetyltransferase(ReactionClass):
    """acetyltransferase activity."""

    GO_ID = "GO:0016407"
    CHILD_CLASSES = (
        NTerminalAminoAcidNAlphaAcetyltransferaseNatA,
        NTerminalMethionineNAlphaAcetyltransferaseNatB,
        NTerminalMethionineNAlphaAcetyltransferaseNatC,
        NTerminalMethionineNAlphaAcetyltransferaseNatE,
        PolysialicAcidOAcetyltransferase,
    )

    @classmethod
    def supports_evaluation(cls, reaction: Reaction) -> bool:
        return aggregate_explicit_child_supports_evaluation(reaction, cls.CHILD_CLASSES)

    def check_membership_impl(self, reaction: Reaction):
        return aggregate_explicit_child_membership(
            reaction,
            self.CHILD_CLASSES,
            "Acetyltransferase",
        )
