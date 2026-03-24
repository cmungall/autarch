"""Focused tests for batch 11 grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NAD_PLUS,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.hydrolase_hydrolyzing_o_glycosyl_compounds import (
    HydrolaseHydrolyzingOGlycosylCompounds,
)
from autarch.ontology.hydroxymethyl_formyl_and_related_transferase import (
    HydroxymethylFormylAndRelatedTransferase,
)
from autarch.ontology.ligase_forming_carbon_sulfur_bonds import (
    LigaseFormingCarbonSulfurBonds,
)
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_iron_sulfur_protein_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor,
)
from autarch.ontology.phosphotransferase_alcohol_group_as_acceptor import (
    PhosphotransferaseAlcoholGroupAsAcceptor,
)

CHEBI_OXIDIZED_2FE2S = "CHEBI:33737"
CHEBI_REDUCED_2FE2S = "CHEBI:33738"
CHEBI_BENZOATE = "CHEBI:16150"
CHEBI_BENZOYL_COA = "CHEBI:57369"
CHEBI_LACTOSE = "CHEBI:17716"
CHEBI_D_GLUCOSE = "CHEBI:4167"
CHEBI_BETA_D_GALACTOSE = "CHEBI:27667"
CHEBI_ADENOSINE = "CHEBI:16335"
CHEBI_ADENINE = "CHEBI:16708"
CHEBI_D_RIBOSE = "CHEBI:47013"
CHEBI_PYRIDOXAL = "CHEBI:17310"
CHEBI_PYRIDOXAL_PHOSPHATE = "CHEBI:597326"
CHEBI_BICARBONATE = "CHEBI:17544"
CHEBI_AMMONIUM = "CHEBI:28938"
CHEBI_CARBAMOYL_PHOSPHATE = "CHEBI:58228"
CHEBI_METHYLENE_THF = "CHEBI:15636"
CHEBI_THF = "CHEBI:57453"
CHEBI_GLYCINE = "CHEBI:57305"
CHEBI_L_SERINE = "CHEBI:33384"
CHEBI_CATECHOL = "CHEBI:17847"
CHEBI_GUAIACOL = "CHEBI:18199"


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_271_positive_small_molecule_kinase() -> None:
    cls = PhosphotransferaseAlcoholGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_PYRIDOXAL,
                "smiles": "[H]C(=O)C1=C(CO)C=NC(C)=C1O",
                "name": "pyridoxal",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_PYRIDOXAL_PHOSPHATE,
                "smiles": "[H]C(=O)C1=C(COP(=O)([O-])[O-])C=NC(C)=C1O",
                "name": "pyridoxal phosphate",
            },
            {
                "chebi_id": CHEBI_ADP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ADP",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "alcohol" in result.explanation.lower() or "kinase" in result.explanation.lower()


def test_ec_271_rejects_carboxyl_phosphorylation() -> None:
    cls = PhosphotransferaseAlcoholGroupAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_BICARBONATE,
                "smiles": "O=C([O-])O",
                "name": "bicarbonate",
            },
            {
                "chebi_id": CHEBI_AMMONIUM,
                "smiles": "[H][N+]([H])([H])[H]",
                "name": "ammonium",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_CARBAMOYL_PHOSPHATE,
                "smiles": "NC(=O)OP(=O)([O-])[O-]",
                "name": "carbamoyl phosphate",
            },
            {
                "chebi_id": CHEBI_ADP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ADP",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0004553_positive_o_glycoside_hydrolase() -> None:
    cls = HydrolaseHydrolyzingOGlycosylCompounds()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_LACTOSE,
                "smiles": "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O",
                "name": "lactose",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {
                "chebi_id": CHEBI_BETA_D_GALACTOSE,
                "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
                "name": "beta-D-galactose",
            },
            {
                "chebi_id": CHEBI_D_GLUCOSE,
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "D-glucose",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "glycos" in result.explanation.lower()


def test_go_0004553_rejects_n_glycoside_hydrolase() -> None:
    cls = HydrolaseHydrolyzingOGlycosylCompounds()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_ADENOSINE,
                "smiles": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
                "name": "adenosine",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {
                "chebi_id": CHEBI_ADENINE,
                "smiles": "Nc1ncnc2ncnc12",
                "name": "adenine",
            },
            {
                "chebi_id": CHEBI_D_RIBOSE,
                "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H]1O",
                "name": "D-ribose",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0016625_positive_ferredoxin_oxo_redox() -> None:
    cls = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_OXIDIZED_2FE2S,
                "smiles": "S1[Fe+]S[Fe+]1",
                "name": "oxidized [2Fe-2S] ferredoxin",
                "count": 2,
            },
            {
                "chebi_id": CHEBI_OXIDIZED_2FE2S,
                "smiles": "S1[Fe+]S[Fe+]1",
                "name": "oxidized [2Fe-2S] ferredoxin",
                "count": 2,
            },
            {
                "chebi_id": "CHEBI:15361",
                "smiles": "CC(=O)C(=O)[O-]",
                "name": "pyruvate",
            },
            {
                "chebi_id": CHEBI_COA,
                "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS",
                "name": "coenzyme A",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_REDUCED_2FE2S,
                "smiles": "S1[Fe]S[Fe+]1",
                "name": "reduced [2Fe-2S] ferredoxin",
                "count": 2,
            },
            {
                "chebi_id": CHEBI_REDUCED_2FE2S,
                "smiles": "S1[Fe]S[Fe+]1",
                "name": "reduced [2Fe-2S] ferredoxin",
                "count": 2,
            },
            {
                "chebi_id": "CHEBI:57288",
                "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "acetyl-CoA",
            },
            {"chebi_id": "CHEBI:16526", "smiles": "O=C=O", "name": "carbon dioxide"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "iron-sulfur" in result.explanation.lower() or "ferredoxin" in result.explanation.lower()


def test_go_0016625_rejects_nad_dehydrogenase() -> None:
    cls = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsIronSulfurProteinAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "CC=O",
                "name": "acetaldehyde",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {
                "smiles": "CC(=O)[O-]",
                "name": "acetate",
            },
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0016742_positive_folate_one_carbon_transfer() -> None:
    cls = HydroxymethylFormylAndRelatedTransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_METHYLENE_THF,
                "smiles": "[H][C@]12CNC3=C(C(=O)NC(N)=N3)N1CN(C1=CC=C(C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])C=C1)C2",
                "name": "5,10-methylene tetrahydrofolate",
            },
            {
                "chebi_id": CHEBI_GLYCINE,
                "smiles": "[NH3+]CC(=O)[O-]",
                "name": "glycine",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            {
                "chebi_id": CHEBI_THF,
                "smiles": "NC1=NC2=C(N[C@@H](CNC3=CC=C(C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])C=C3)CN2)C(=O)N1",
                "name": "tetrahydrofolate",
            },
            {
                "chebi_id": CHEBI_L_SERINE,
                "smiles": "[NH3+][C@@H](CO)C(=O)[O-]",
                "name": "L-serine",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "folate" in result.explanation.lower() or "one-carbon" in result.explanation.lower()


def test_go_0016742_rejects_sam_methyltransferase() -> None:
    cls = HydroxymethylFormylAndRelatedTransferase()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_SAM,
                "name": "SAM",
            },
            {
                "chebi_id": CHEBI_CATECHOL,
                "smiles": "Oc1ccccc1O",
                "name": "catechol",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_SAH,
                "name": "SAH",
            },
            {
                "chebi_id": CHEBI_GUAIACOL,
                "smiles": "COc1ccccc1O",
                "name": "guaiacol",
            },
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_go_0016877_positive_acyl_thioester_ligase() -> None:
    cls = LigaseFormingCarbonSulfurBonds()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_BENZOATE,
                "smiles": "O=C([O-])C1=CC=CC=C1",
                "name": "benzoate",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
            {
                "chebi_id": CHEBI_COA,
                "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS",
                "name": "coenzyme A",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_BENZOYL_COA,
                "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1=CC=CC=C1",
                "name": "benzoyl-CoA",
            },
            {
                "chebi_id": CHEBI_AMP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "AMP",
            },
            {
                "chebi_id": CHEBI_DIPHOSPHATE,
                "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                "name": "diphosphate",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "carbon-sulfur" in result.explanation.lower() or "thioester" in result.explanation.lower()


def test_go_0016877_rejects_amide_ligase() -> None:
    cls = LigaseFormingCarbonSulfurBonds()
    reaction = _reaction(
        left=[
            {
                "chebi_id": CHEBI_BENZOATE,
                "smiles": "O=C([O-])C1=CC=CC=C1",
                "name": "benzoate",
            },
            {
                "chebi_id": CHEBI_ATP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "ATP",
            },
            {
                "chebi_id": CHEBI_AMMONIUM,
                "smiles": "[H][N+]([H])([H])[H]",
                "name": "ammonium",
            },
        ],
        right=[
            {
                "smiles": "NC(=O)C1=CC=CC=C1",
                "name": "benzamide",
            },
            {
                "chebi_id": CHEBI_AMP,
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                "name": "AMP",
            },
            {
                "chebi_id": CHEBI_DIPHOSPHATE,
                "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                "name": "diphosphate",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
