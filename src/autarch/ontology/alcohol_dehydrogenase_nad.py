"""alcohol dehydrogenase NAD."""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import Moiety, is_aldehyde, is_carboxylic_acid
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.pattern_dsl import var, match_patterns
from autarch.molecules import (
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    h_plus,
    nad_plus,
    nadh,
    p,
)


class AlcoholDehydrogenaseNAD(OxidoreductaseActingOnTheCHOHGroupOfDonorsNADOrNADPAsAcceptor):
    """alcohol dehydrogenase NAD

    Examples:
    - Ethanol dehydrogenase: ethanol + NAD+ → acetaldehyde + NADH + H+
    - Methanol dehydrogenase: methanol + NAD+ → formaldehyde + NADH + H+
    - Lactate dehydrogenase: lactate + NAD+ → pyruvate + NADH + H+

    Scope note:
    GO:0004022 covers alcohol <-> aldehyde/ketone interconversion with NAD+/NADH.
    This classifier focuses on that canonical transformation.
    """

    GO_ID = "GO:0004022"  # alcohol dehydrogenase (NAD+) activity
    EC_NUMBER_PREFIX = "1.1.1.1"  # alcohol dehydrogenase (NAD+) activity
    KETONE_PATTERN = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")

    # Define patterns using DSL
    PATTERNS: list[Reaction] = [
        # Primary alcohol → aldehyde: R-CH2OH + NAD+ → R-CHO + NADH + H+
        var("alcohol") + p(nad_plus) >> var("aldehyde") + p(nadh) + p(h_plus),
        # Secondary alcohol → ketone: R-CHOH-R' + NAD+ → R-CO-R' + NADH + H+
        var("secondary_alcohol") + p(nad_plus) >> var("ketone") + p(nadh) + p(h_plus),
        # Reverse (reduction): aldehyde/ketone + NADH + H+ → alcohol + NAD+
        var("carbonyl") + p(nadh) + p(h_plus) >> var("alcohol") + p(nad_plus),
    ]

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an alcohol dehydrogenase using pattern matching.

        Strategy:
        1. Must be a dehydrogenase (parent class)
        2. Must involve alcohol oxidation or aldehyde/ketone reduction
        3. Must use NAD+/NADH system
        4. Look for characteristic alcohol/carbonyl substrates
        """
        # First check if it's a dehydrogenase
        parent_result = super().check_membership_impl(reaction)
        if not parent_result.is_member and "No clear NAD(P)+/NAD(P)H conversion pattern" not in parent_result.explanation:
            return ClassificationResult(
                is_member=False,
                explanation=f"Not a dehydrogenase: {parent_result.explanation}",
            )

        # GO:0004022 is NAD+ specific; reject NADP+ dependent variants.
        has_nadp = any(
            participant.chebi_id in {CHEBI_NADP_PLUS, CHEBI_NADPH}
            for participant in reaction.left_participants + reaction.right_participants
        )
        if has_nadp:
            return ClassificationResult(
                is_member=False,
                explanation="Uses NADP(H) cofactor - not NAD+-specific alcohol dehydrogenase",
            )

        has_nad_plus_left = any(
            p.chebi_id == CHEBI_NAD_PLUS for p in reaction.left_participants
        )
        has_nad_plus_right = any(
            p.chebi_id == CHEBI_NAD_PLUS for p in reaction.right_participants
        )
        has_nadh_left = any(p.chebi_id == CHEBI_NADH for p in reaction.left_participants)
        has_nadh_right = any(
            p.chebi_id == CHEBI_NADH for p in reaction.right_participants
        )
        is_forward = has_nad_plus_left and has_nadh_right
        is_reverse = has_nadh_left and has_nad_plus_right
        if not (is_forward or is_reverse):
            return ClassificationResult(
                is_member=False,
                explanation="No NAD+/NADH conversion direction for alcohol dehydrogenase",
            )

        ignored_ids = {
            CHEBI_NAD_PLUS,
            CHEBI_NADH,
            CHEBI_NADP_PLUS,
            CHEBI_NADPH,
            CHEBI_H_PLUS,
        }
        left_non_cofactors = [
            p for p in reaction.left_participants if p.chebi_id not in ignored_ids
        ]
        right_non_cofactors = [
            p for p in reaction.right_participants if p.chebi_id not in ignored_ids
        ]

        # Alcohol dehydrogenase (GO:0004022) is a simple 1->1 alcohol/carbonyl conversion.
        if len(left_non_cofactors) != 1 or len(right_non_cofactors) != 1:
            return ClassificationResult(
                is_member=False,
                explanation="Not a simple 1->1 alcohol dehydrogenase substrate/product pair",
            )

        noncofactors = left_non_cofactors + right_non_cofactors
        has_phosphorylated = any(p.is_phosphorylated() for p in noncofactors)
        has_carboxylate = any(
            p.smiles and is_carboxylic_acid(p.smiles, p.chebi_id) for p in noncofactors
        )
        has_thioester = any(p.is_thioester() for p in noncofactors)
        if has_phosphorylated or has_carboxylate or has_thioester:
            return ClassificationResult(
                is_member=False,
                explanation="Phosphorylated/carboxylate/thioester substrate class - not GO:0004022 alcohol dehydrogenase",
            )

        has_alcohol_left = any(self._is_alcohol_like(p) for p in left_non_cofactors)
        has_alcohol_right = any(self._is_alcohol_like(p) for p in right_non_cofactors)
        has_aldehyde_left = any(self._is_aldehyde_like(p) for p in left_non_cofactors)
        has_aldehyde_right = any(self._is_aldehyde_like(p) for p in right_non_cofactors)
        has_ketone_left = any(self._is_ketone_like(p) for p in left_non_cofactors)
        has_ketone_right = any(self._is_ketone_like(p) for p in right_non_cofactors)

        is_alcohol_oxidation = is_forward and has_alcohol_left and (
            has_aldehyde_right or has_ketone_right
        )
        is_carbonyl_reduction = is_reverse and has_alcohol_right and (
            has_aldehyde_left or has_ketone_left
        )

        if not (is_alcohol_oxidation or is_carbonyl_reduction):
            return ClassificationResult(
                is_member=False,
                explanation="No alcohol ⇌ carbonyl conversion detected",
            )

        # Try to match patterns
        match = match_patterns(reaction, self.PATTERNS, strict=False)

        explanation = "Alcohol dehydrogenase: alcohol ⇌ carbonyl"

        # Check direction
        if is_alcohol_oxidation:
            explanation += " [oxidation]"
        elif is_carbonyl_reduction:
            explanation += " [reduction]"

        explanation += " [NAD-dependent]"

        if match.matched:
            explanation += " [pattern matched]"

        return ClassificationResult(is_member=True, explanation=explanation)

    @classmethod
    def _is_alcohol_like(cls, participant: Participant) -> bool:
        return participant.has_moiety(Moiety.HYDROXYL)

    @classmethod
    def _is_aldehyde_like(cls, participant: Participant) -> bool:
        return bool(
            participant.smiles and is_aldehyde(participant.smiles, participant.chebi_id)
        )

    @classmethod
    def _is_ketone_like(cls, participant: Participant) -> bool:
        if not participant.smiles:
            return False
        mol = participant.get_mol()
        if mol is None:
            return False
        if cls.KETONE_PATTERN is not None and mol.HasSubstructMatch(cls.KETONE_PATTERN):
            return True
        return False
