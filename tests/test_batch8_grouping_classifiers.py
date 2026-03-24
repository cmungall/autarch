"""Focused tests for batch 8 grouping classifiers."""

from __future__ import annotations

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CO2,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.carbon_nitrogen_ligase_with_glutamine_as_amido_n_donor import (
    CarbonNitrogenLigaseWithGlutamineAsAmidoNDonor,
)
from autarch.ontology.ether_hydrolase import EtherHydrolase
from autarch.ontology.hydrolase_acting_on_halide_bonds_in_c_halide_compounds import (
    HydrolaseActingOnHalideBondsInCHalideCompounds,
)
from autarch.ontology.oxidoreductase_acting_on_a_sulfur_group_of_donors_disulfide_as_acceptor import (
    OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_single_donors_with_incorporation_of_molecular_oxygen_incorporation_of_one_atom_of_oxygen_internal_monooxygenases_or_internal_mixed_function_oxidases import (
    OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfOneAtomOfOxygenInternalMonooxygenasesOrInternalMixedFunctionOxidases,
)


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    inchi: str
    name: str
    count: int
    polymer_type: PolymerType


ATP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])"
    "[C@@H](O)[C@H]1O"
)
ADP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
)
PHOSPHATE_SMILES = "O=P([O-])([O-])[O-]"


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_184_positive_disulfide_acceptor_redox() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC(CCSSCCC(N)C(=O)O)C(=O)O",
                "name": "L-homocystine",
            },
            {
                "smiles": "N[C@@H](CCC(=O)NCCS)C(=O)O",
                "name": "glutathione",
                "count": 2,
            },
        ],
        right=[
            {
                "smiles": "N[C@@H](CCS)C(=O)O",
                "name": "L-homocysteine",
                "count": 2,
            },
            {
                "smiles": "N[C@@H](CCC(=O)NCCSSC[C@H](N)C(=O)O)C(=O)O",
                "name": "glutathione disulfide",
            },
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "disulfide" in result.explanation.lower()


def test_ec_184_rejects_nadph_linked_sulfur_redox() -> None:
    cls = OxidoreductaseActingOnASulfurGroupOfDonorsDisulfideAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "N[C@@H](CCC(=O)NCCSSC[C@H](N)C(=O)O)C(=O)O",
                "name": "glutathione disulfide",
            },
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "smiles": "N[C@@H](CCC(=O)NCCS)C(=O)O",
                "name": "glutathione",
                "count": 2,
            },
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_332_positive_epoxide_hydrolysis() -> None:
    cls = EtherHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "C1OC1c1ccccc1", "name": "styrene oxide"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "OC(CO)c1ccccc1", "name": "phenylethane-1,2-diol"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "ether" in result.explanation.lower()


def test_ec_332_rejects_simple_ester_hydrolysis() -> None:
    cls = EtherHydrolase()
    reaction = _reaction(
        left=[
            {"smiles": "CC(=O)OC", "name": "methyl acetate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"smiles": "CO", "name": "methanol"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_635_positive_glutamine_dependent_amidation() -> None:
    cls = CarbonNitrogenLigaseWithGlutamineAsAmidoNDonor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:78516",
                "name": "L-aspartyl-tRNA(Asn)",
                "polymer_type": PolymerType.TRNA,
                "smiles": "*P(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OC(=O)[C@@H]([NH3+])CC(=O)[O-]",
            },
            {
                "chebi_id": "CHEBI:58359",
                "smiles": "NC(=O)CC[C@H]([NH3+])C(=O)[O-]",
                "name": "L-glutamine",
            },
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {
                "chebi_id": "CHEBI:78515",
                "name": "L-asparaginyl-tRNA(Asn)",
                "polymer_type": PolymerType.TRNA,
                "smiles": "*P(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OC(=O)[C@@H](N)CC(N)=O",
            },
            {
                "chebi_id": "CHEBI:29985",
                "smiles": "[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]",
                "name": "L-glutamate",
            },
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "glutamine" in result.explanation.lower()


def test_ec_635_rejects_ammonia_ligase_branch() -> None:
    cls = CarbonNitrogenLigaseWithGlutamineAsAmidoNDonor()
    reaction = _reaction(
        left=[
            {"smiles": "[NH3+]CCC[C@H]([NH3+])C(=O)[O-]", "name": "L-lysine"},
            {"smiles": "O=C([O-])CC[C@H]([NH3+])C(=O)[O-]", "name": "L-glutamate"},
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
            {"smiles": "N", "name": "ammonia"},
        ],
        right=[
            {
                "smiles": "[NH3+]CCC[C@H](NC(=O)CC[C@H]([NH3+])C(=O)[O-])C(=O)[O-]",
                "name": "glutamyl-lysine",
            },
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
            {"chebi_id": CHEBI_PHOSPHATE, "smiles": PHOSPHATE_SMILES, "name": "phosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_11312_positive_internal_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfOneAtomOfOxygenInternalMonooxygenasesOrInternalMixedFunctionOxidases()
    )
    reaction = _reaction(
        left=[
            {
                "smiles": "[NH3+][C@@H](CC1=CC=CC=C1)C(=O)[O-]",
                "name": "L-phenylalanine",
            },
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
        ],
        right=[
            {"smiles": "NC(=O)CC1=CC=CC=C1", "name": "2-phenylacetamide"},
            {"chebi_id": CHEBI_CO2, "smiles": "O=C=O", "name": "carbon dioxide"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "internal monooxygenase" in result.explanation.lower()


def test_ec_11312_rejects_nadph_monooxygenase() -> None:
    cls = (
        OxidoreductaseActingOnSingleDonorsWithIncorporationOfMolecularOxygenIncorporationOfOneAtomOfOxygenInternalMonooxygenasesOrInternalMixedFunctionOxidases()
    )
    reaction = _reaction(
        left=[
            {"smiles": "c1ccc(cc1)CC(=O)[O-]", "name": "phenylacetate"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "c1cc(O)c(cc1)CC(=O)[O-]", "name": "hydroxyphenylacetate"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP(+)"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_381_positive_haloalkane_hydrolysis() -> None:
    cls = HydrolaseActingOnHalideBondsInCHalideCompounds()
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCl", "name": "1-chlorohexane"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CCCCCCO", "name": "hexan-1-ol"},
            {"chebi_id": "CHEBI:17996", "smiles": "[Cl-]", "name": "chloride"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "halide" in result.explanation.lower()


def test_ec_381_rejects_reductive_dehalogenation() -> None:
    cls = HydrolaseActingOnHalideBondsInCHalideCompounds()
    reaction = _reaction(
        left=[
            {"smiles": "CCCCCCCl", "name": "1-chlorohexane"},
            {"smiles": "[H][H]", "name": "dihydrogen"},
        ],
        right=[
            {"smiles": "CCCCCC", "name": "hexane"},
            {"chebi_id": "CHEBI:17996", "smiles": "[Cl-]", "name": "chloride"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False
