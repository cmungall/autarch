"""carbon-carbon lyase.

Catalysis of the cleavage of C-C bonds by other means than by hydrolysis or
oxidation, or conversely adding a group to a double bond.
"""

from autarch.ontology.carbon_carbon_lyases import CarbonCarbonLyases


class CarbonCarbonLyase(CarbonCarbonLyases):
    """carbon-carbon lyase.

    Catalysis of the cleavage of C-C bonds by other means than by hydrolysis or
    oxidation, or conversely adding a group to a double bond.
    """

    GO_ID = "GO:0016830"
    EC_NUMBER_PREFIX = "4.1.-.-"
