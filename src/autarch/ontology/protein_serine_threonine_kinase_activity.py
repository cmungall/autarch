"""protein serine/threonine kinase activity.

Catalysis of the reactions: ATP + protein serine = ADP + protein serine
phosphate, and ATP + protein threonine = ADP + protein threonine phosphate.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.protein_kinase_activity import ProteinKinaseActivity


class ProteinSerineThreonineKinaseActivity(ProteinKinaseActivity):
    """protein serine/threonine kinase activity.

    Catalysis of the reactions: ATP + protein serine = ADP + protein serine
    phosphate, and ATP + protein threonine = ADP + protein threonine
    phosphate.
    """

    GO_ID = "GO:0004674"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        parent = super().check_membership_impl(reaction)
        if not parent.is_member:
            return parent
        return ClassificationResult(
            is_member=True,
            explanation="Protein serine/threonine kinase activity: ATP-dependent phosphorylation of a protein substrate in the serine/threonine kinase branch",
        )
