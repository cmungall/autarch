"""ligase activity, forming carbon-sulfur bonds.

Catalysis of the joining of two molecules via a carbon-sulfur bond, with the
concomitant hydrolysis of the diphosphate bond in ATP or a similar triphosphate.
"""

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.ligase import Ligase

CHEBI_ALT_AMP = "CHEBI:456215"
CHEBI_ALT_ADP = "CHEBI:456216"
CHEBI_ALT_ATP = "CHEBI:30616"
CHEBI_ALT_COA = "CHEBI:57287"
CHEBI_ALT_DIPHOSPHATE = "CHEBI:33019"
CHEBI_ALT_GDP = "CHEBI:58189"
CHEBI_ALT_GTP = "CHEBI:37565"
CHEBI_ALT_PHOSPHATE = "CHEBI:16838"


class LigaseFormingCarbonSulfurBonds(Ligase):
    """ligase activity, forming carbon-sulfur bonds.

    Catalysis of the joining of two molecules via a carbon-sulfur bond, with
    the concomitant hydrolysis of the diphosphate bond in ATP or a similar
    triphosphate.
    """

    GO_ID = "GO:0016877"
    EC_NUMBER_PREFIX = "6.2.-.-"
    ATP_CHEBIS = {CHEBI_ATP, CHEBI_ALT_ATP}
    AMP_CHEBIS = {CHEBI_AMP, CHEBI_ALT_AMP}
    ADP_CHEBIS = {CHEBI_ALT_ADP}
    GTP_CHEBIS = {CHEBI_GTP, CHEBI_ALT_GTP}
    GDP_CHEBIS = {CHEBI_GDP, CHEBI_ALT_GDP}
    PHOSPHATE_CHEBIS = {CHEBI_PHOSPHATE, CHEBI_ALT_PHOSPHATE}
    DIPHOSPHATE_CHEBIS = {CHEBI_DIPHOSPHATE, CHEBI_ALT_DIPHOSPHATE}
    COA_CHEBIS = {CHEBI_COA, CHEBI_ALT_COA}
    CARRIER_PROTEIN_CHEBIS = {
        "CHEBI:64479",
        "CHEBI:133243",
        "CHEBI:138622",
        "CHEBI:29950",
        "CHEBI:133479",
        "CHEBI:82683",
        "CHEBI:137976",
        "CHEBI:78446",
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for ATP- or GTP-dependent formation of a carbon-sulfur bond."""
        if any(participant.chebi_id == CHEBI_H2O for participant in reaction.left_participants):
            return ClassificationResult(
                is_member=False,
                explanation="Uses water as a substrate - not carbon-sulfur ligase chemistry",
            )

        has_atp = any(participant.chebi_id in self.ATP_CHEBIS for participant in reaction.left_participants)
        has_gtp = any(participant.chebi_id in self.GTP_CHEBIS for participant in reaction.left_participants)
        has_amp = any(participant.chebi_id in self.AMP_CHEBIS for participant in reaction.right_participants)
        has_adp = any(participant.chebi_id in self.ADP_CHEBIS for participant in reaction.right_participants)
        has_gdp = any(participant.chebi_id in self.GDP_CHEBIS for participant in reaction.right_participants)
        has_phosphate = any(
            participant.chebi_id in self.PHOSPHATE_CHEBIS for participant in reaction.right_participants
        )
        has_diphosphate = any(
            participant.chebi_id in self.DIPHOSPHATE_CHEBIS for participant in reaction.right_participants
        )

        has_nucleotide_coupling = (has_atp and ((has_amp and has_diphosphate) or (has_adp and has_phosphate))) or (
            has_gtp and has_gdp and has_phosphate
        )
        if not has_nucleotide_coupling:
            return ClassificationResult(
                is_member=False,
                explanation="No ATP/GTP coupling pattern consistent with carbon-sulfur ligase chemistry",
            )

        has_acid_substrate = any(
            participant.has_moiety(Moiety.CARBOXYL) and not participant.is_thioester()
            for participant in reaction.left_participants
        )
        if not has_acid_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No acid substrate detected",
            )

        has_thiol_substrate = any(self._is_thiol_substrate(participant) for participant in reaction.left_participants)
        if not has_thiol_substrate:
            return ClassificationResult(
                is_member=False,
                explanation="No thiol substrate detected",
            )

        has_thioester_product = any(participant.is_thioester() for participant in reaction.right_participants)
        if not has_thioester_product:
            return ClassificationResult(
                is_member=False,
                explanation="No thioester product detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Carbon-sulfur ligase: ATP-dependent formation of a thioester or related carbon-sulfur bond",
        )

    @classmethod
    def _is_thiol_substrate(cls, participant: Participant) -> bool:
        if participant.chebi_id in cls.COA_CHEBIS | cls.CARRIER_PROTEIN_CHEBIS:
            return True
        if participant.is_thioester():
            return False
        mol = participant.get_mol()
        if mol is None:
            return False
        return any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
