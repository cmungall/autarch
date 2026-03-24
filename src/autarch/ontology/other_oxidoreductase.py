"""Other oxidoreductase reaction classification.

Catch-all for oxidoreductases (EC 1.97) that do not fit into more specific
subcategories. These are redox reactions with unusual or poorly characterized
electron donors/acceptors.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase


class OtherOxidoreductase(Oxidoreductase):
    """Catalysis of an oxidation-reduction (redox) reaction, a reversible chemical reaction in which the oxidation state of an atom or atoms within a molecule is altered."""

    GO_ID = "GO:0016491"  # oxidoreductase activity (same as parent - catch-all)
    EC_NUMBER_PREFIX = "1.97.-.-"  # Other oxidoreductases

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an oxidoreductase.

        This is a catch-all class (EC 1.97) for oxidoreductases that don't fit
        more specific subcategories. It simply delegates to the parent
        Oxidoreductase classifier.
        """
        return super().check_membership_impl(reaction)
