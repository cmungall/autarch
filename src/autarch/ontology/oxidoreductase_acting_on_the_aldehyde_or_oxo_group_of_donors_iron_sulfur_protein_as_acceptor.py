"""oxidoreductase activity, acting on the aldehyde or oxo group of donors, iron-sulfur protein as acceptor.

Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde or
ketone (oxo) group acts as a hydrogen or electron donor and reduces an
iron-sulfur protein.
"""

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, Reaction
from autarch.moiety import is_aldehyde
from autarch.molecules import (
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
)
from autarch.ontology.oxidoreductase import Oxidoreductase

CHEBI_CO = "CHEBI:17245"


class OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor(
    Oxidoreductase
):
    """oxidoreductase activity, acting on the aldehyde or oxo group of donors, iron-sulfur protein as acceptor.

    Catalysis of an oxidation-reduction (redox) reaction in which an aldehyde
    or ketone (oxo) group acts as a hydrogen or electron donor and reduces an
    iron-sulfur protein.
    """

    GO_ID = "GO:0016625"
    EC_NUMBER_PREFIX = "1.2.7.-"
    OXIDIZED_IRON_SULFUR_CHEBIS = {"CHEBI:33737", "CHEBI:33722"}
    REDUCED_IRON_SULFUR_CHEBIS = {"CHEBI:33738", "CHEBI:33723"}
    KETO_ACID_PATTERN = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)")

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for iron-sulfur-dependent oxidation of aldehyde/oxo donors."""
        has_oxidized_acceptor = any(
            participant.chebi_id in self.OXIDIZED_IRON_SULFUR_CHEBIS
            for participant in reaction.left_participants
        )
        has_reduced_product = any(
            participant.chebi_id in self.REDUCED_IRON_SULFUR_CHEBIS
            for participant in reaction.right_participants
        )
        if not (has_oxidized_acceptor and has_reduced_product):
            return ClassificationResult(
                is_member=False,
                explanation="No oxidized iron-sulfur acceptor -> reduced iron-sulfur product signature",
            )

        if any(
            participant.chebi_id in {
                CHEBI_NAD_PLUS,
                CHEBI_NADP_PLUS,
                CHEBI_NADH,
                CHEBI_NADPH,
                CHEBI_O2,
                CHEBI_H2O2,
            }
            for participant in reaction.all_participants()
        ):
            return ClassificationResult(
                is_member=False,
                explanation="Competing nicotinamide or oxygen acceptor chemistry indicates a different oxidoreductase branch",
            )

        left_core = [
            participant
            for participant in reaction.left_participants
            if participant.chebi_id not in self.OXIDIZED_IRON_SULFUR_CHEBIS | {CHEBI_H2O, CHEBI_H_PLUS}
        ]
        right_core = [
            participant
            for participant in reaction.right_participants
            if participant.chebi_id not in self.REDUCED_IRON_SULFUR_CHEBIS | {CHEBI_H_PLUS}
        ]
        if not left_core or not right_core:
            return ClassificationResult(
                is_member=False,
                explanation="No non-cofactor donor/product pair to evaluate",
            )

        donor_kinds = {kind for participant in left_core for kind in self._donor_kinds(participant)}
        if not donor_kinds:
            return ClassificationResult(
                is_member=False,
                explanation="No aldehyde or oxo donor detected on the substrate side",
            )

        has_thioester_product = any(participant.is_thioester() for participant in right_core)
        has_carboxylate_product = any(participant.is_carboxylic_acid() for participant in right_core)
        has_co2_product = any(participant.chebi_id == CHEBI_CO2 for participant in right_core)
        has_oxidized_product = has_thioester_product or (
            {"aldehyde", "formate", "co"} & donor_kinds and (has_carboxylate_product or has_co2_product)
        ) or ("keto_acid" in donor_kinds and (has_thioester_product or has_co2_product))
        if not has_oxidized_product:
            return ClassificationResult(
                is_member=False,
                explanation="No product consistent with aldehyde/oxo oxidation was detected",
            )

        return ClassificationResult(
            is_member=True,
            explanation="Iron-sulfur oxidoreductase: aldehyde/oxo donor coupled to reduction of an iron-sulfur acceptor",
        )

    @classmethod
    def _donor_kinds(cls, participant: Participant) -> set[str]:
        kinds: set[str] = set()
        if participant.chebi_id == CHEBI_CO:
            kinds.add("co")
            return kinds
        if not participant.smiles:
            return kinds
        if is_aldehyde(participant.smiles, participant.chebi_id):
            kinds.add("aldehyde")
        mol = participant.get_mol()
        if mol is not None and cls.KETO_ACID_PATTERN is not None and mol.HasSubstructMatch(cls.KETO_ACID_PATTERN):
            kinds.add("keto_acid")
        if mol is not None and cls._is_formyl_heteroatom_adduct(mol):
            kinds.add("formate")
        return kinds

    @staticmethod
    def _is_formyl_heteroatom_adduct(mol: Chem.Mol) -> bool:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            oxygen_neighbors = 0
            hetero_single = False
            hydrogen_count = atom.GetTotalNumHs()
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxygen_neighbors += 1
                if neighbor.GetAtomicNum() in {7, 8, 16} and bond.GetBondType() == Chem.BondType.SINGLE:
                    hetero_single = True
            if oxygen_neighbors == 1 and hetero_single and hydrogen_count >= 1:
                return True
        return False
