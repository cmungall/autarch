"""acting on NADH or NADPH."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.nadh_or_nadph_dehydrogenase_quinone import NADHOrNADPHDehydrogenaseQuinone
from autarch.ontology.oxidoreductase_acting_on_nadh_or_nadph_quinone_or_similar_compound_as_acceptor import OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor


class ActingOnNADHOrNADPH(ExplicitEcAggregate):
    """acting on NADH or NADPH."""

    EC_NUMBER_PREFIX = '1.6.-.-'
    CHILD_CLASSES = (
        NADHOrNADPHDehydrogenaseQuinone,
        OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor,
    )
