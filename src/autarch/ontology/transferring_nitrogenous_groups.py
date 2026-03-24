"""transferring nitrogenous groups."""

from autarch.ontology.ec_prefix_aggregate import ExplicitEcAggregate
from autarch.ontology.aromatic_amino_acid_transaminase import AromaticAminoAcidTransaminase
from autarch.ontology.transaminase import Transaminase


class TransferringNitrogenousGroups(ExplicitEcAggregate):
    """transferring nitrogenous groups."""

    EC_NUMBER_PREFIX = '2.6.-.-'
    CHILD_CLASSES = (
        AromaticAminoAcidTransaminase,
        Transaminase,
    )
