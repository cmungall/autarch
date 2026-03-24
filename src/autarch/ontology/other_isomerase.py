"""Other isomerase reaction classification.

EC 5.99: Other isomerases that do not fit into the standard isomerase
subcategories. This is a catch-all for unusual isomerization reactions.
"""

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.isomerase import Isomerase


class OtherIsomerase(Isomerase):
    """other isomerase

    EC 5.99.x.x includes miscellaneous isomerases that do not fit
    standard subcategories.

    Examples:
    - Thiocyanate isomerase
    - DNA adenine methyltransferase (isomerization component)
    """

    GO_ID = "GO:0016853"  # isomerase activity (same as parent - catch-all)
    EC_NUMBER_PREFIX = "5.99.-.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an 'other' isomerase (EC 5.99).

        Delegates to parent Isomerase for general isomerase classification.
        """
        return super().check_membership_impl(reaction)
