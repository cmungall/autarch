"""transferase transferring sulfur-containing groups.

Catalysis of the transfer of a sulfur-containing group from one compound
(donor) to another (acceptor).
"""

from autarch.ontology.transferring_sulfur_containing_groups import (
    TransferringSulfurContainingGroups,
)


class TransferaseTransferringSulfurContainingGroups(
    TransferringSulfurContainingGroups
):
    """transferase transferring sulfur-containing groups.

    Catalysis of the transfer of a sulfur-containing group from one compound
    (donor) to another (acceptor).
    """

    GO_ID = "GO:0016782"
    EC_NUMBER_PREFIX = "2.8.-.-"
