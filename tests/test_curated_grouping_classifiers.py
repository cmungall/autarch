"""Focused tests for chemistry-driven grouping classifiers."""

from typing_extensions import TypedDict

from autarch.datamodel import Participant, Reaction
from autarch.molecules import (
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_CO2,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
)
from autarch.ontology.acid_thiol_ligase import AcidThiolLigase
from autarch.ontology.carboxylic_ester_hydrolase import CarboxylicEsterHydrolase
from autarch.ontology.pyrophosphatase import Pyrophosphatase
from autarch.ontology.oxidoreductase_acting_on_the_aldehyde_or_oxo_group_of_donors_nad_or_nadp_as_acceptor import (
    OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor,
)
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor,
)


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in grouping-classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int


def make_reaction(
    left_participants: list[ParticipantInput],
    right_participants: list[ParticipantInput],
) -> Reaction:
    """Build a reaction for classifier tests."""
    return Reaction(
        left_participants=[Participant(**participant) for participant in left_participants],
        right_participants=[Participant(**participant) for participant in right_participants],
    )


class TestCHOHOxygenAcceptor:
    """Tests for EC 1.1.3 grouping logic."""

    def test_detects_glucose_oxidase_pattern(self):
        classifier = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
                {
                    "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                    "name": "glucose",
                },
            ],
            right_participants=[
                {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
                {
                    "smiles": "OC[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1O",
                    "name": "gluconolactone",
                },
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True
        assert "oxygen acceptor" in result.explanation.lower()

    def test_detects_secondary_alcohol_oxidation(self):
        classifier = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
                {"smiles": "CC(O)C", "name": "isopropanol"},
            ],
            right_participants=[
                {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
                {"smiles": "CC(=O)C", "name": "acetone"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_rejects_amine_oxidation(self):
        classifier = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"chebi_id": CHEBI_O2, "smiles": "O=O", "name": "dioxygen"},
                {"smiles": "CN", "name": "methylamine"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"chebi_id": CHEBI_H2O2, "smiles": "OO", "name": "hydrogen peroxide"},
                {"smiles": "C=O", "name": "formaldehyde"},
                {"smiles": "N", "name": "ammonia"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False


class TestAldehydeOrOxoNADPAcceptor:
    """Tests for EC 1.2 grouping logic."""

    def test_detects_aldehyde_oxidation_to_carboxylate(self):
        classifier = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "[H]C(=O)C1=CC=CC=C1", "name": "benzaldehyde"},
                {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "O=C([O-])C1=CC=CC=C1", "name": "benzoate"},
                {"chebi_id": CHEBI_NADH, "name": "NADH"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True
        assert "aldehyde/oxo" in result.explanation.lower()

    def test_detects_aldehyde_oxidation_to_acyl_phosphate(self):
        classifier = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC=O", "name": "acetaldehyde"},
                {"smiles": "OP(=O)([O-])[O-]", "name": "phosphate"},
                {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
            ],
            right_participants=[
                {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
                {"chebi_id": CHEBI_NADH, "name": "NADH"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_oxo_oxidative_decarboxylation_to_thioester(self):
        classifier = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC(=O)C(=O)[O-]", "name": "pyruvate"},
                {"smiles": "NCCS", "name": "thiol acceptor"},
                {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
            ],
            right_participants=[
                {"smiles": "CC(=O)SCCN", "name": "acetyl thioester"},
                {"chebi_id": CHEBI_CO2, "smiles": "O=C=O", "name": "carbon dioxide"},
                {"chebi_id": CHEBI_NADH, "name": "NADH"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_rejects_alcohol_dehydrogenase_pattern(self):
        classifier = OxidoreductaseActingOnTheAldehydeOrOxoGroupOfDonorsNADOrNADPAsAcceptor()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CCO", "name": "ethanol"},
                {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD(+)"},
            ],
            right_participants=[
                {"smiles": "CC=O", "name": "acetaldehyde"},
                {"chebi_id": CHEBI_NADH, "name": "NADH"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False


class TestPyrophosphatase:
    """Tests for EC 3.6.1 grouping logic."""

    def test_detects_nucleoside_triphosphate_hydrolysis(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "O=C1C=CN([C@H]2C[C@H](O)[C@@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O2)C(=O)N1",
                    "name": "dUTP",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "O=C1C=CN([C@H]2C[C@H](O)[C@@H](COP(=O)([O-])[O-])O2)C(=O)N1",
                    "name": "dUMP",
                },
                {
                    "chebi_id": CHEBI_DIPHOSPHATE,
                    "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                    "name": "diphosphate",
                },
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_adp_ribose_hydrolysis(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "NC1=C2/N=C\\N([C@@H]3O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]4OC(O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H]3O)C2=NC=N1",
                    "name": "ADP-ribose",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "O=P([O-])([O-])OC[C@H]1OC(O)[C@H](O)[C@@H]1O",
                    "name": "ribose 5-phosphate",
                },
                {
                    "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                    "name": "AMP",
                },
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_inorganic_diphosphate_hydrolysis(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "chebi_id": CHEBI_DIPHOSPHATE,
                    "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                    "name": "diphosphate",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_rejects_prenyl_diphosphate_hydrolysis(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "CC(C)=CCC/C(C)=C/COP(=O)([O-])OP(=O)([O-])[O-]",
                    "name": "geranyl diphosphate",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "C=C[C@](C)(O)CCC=C(C)C", "name": "linalool"},
                {
                    "chebi_id": CHEBI_DIPHOSPHATE,
                    "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                    "name": "diphosphate",
                },
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False


class TestCarboxylicEsterHydrolase:
    """Tests for EC 3.1.1 grouping logic."""

    def test_detects_simple_carboxylic_ester_hydrolysis(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "COC(=O)c1ccccc1O", "name": "methyl salicylate"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "CO", "name": "methanol"},
                {"smiles": "O=C([O-])c1ccccc1O", "name": "salicylate"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_reverse_orientation_of_ester_hydrolysis(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CO", "name": "methanol"},
                {"smiles": "O=C([O-])c1ccccc1O", "name": "salicylate"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
            right_participants=[
                {"smiles": "COC(=O)c1ccccc1O", "name": "methyl salicylate"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_lactone_ring_opening(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "O=C1OCCCc2ccccc21", "name": "3,4-dihydrocoumarin"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "O=C([O-])CCCc1ccccc1O", "name": "3-(2-hydroxyphenyl)propanoate"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_tautomerized_lactone_hydrolysis(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "O=C([O-])CC1=CCC(=O)O1", "name": "5-oxo-4,5-dihydro-2-furylacetate"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "O=C([O-])CCC(=O)CC(=O)[O-]", "name": "3-oxoadipate"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_rejects_amide_hydrolysis(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC(=O)N", "name": "acetamide"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "CC(=O)[O-]", "name": "acetate"},
                {"smiles": "N", "name": "ammonia"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False


class TestAcidThiolLigase:
    """Tests for EC 6.2.1 grouping logic."""

    def test_detects_amp_diphosphate_coupling(self):
        classifier = AcidThiolLigase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC(=O)[O-]", "name": "acetate"},
                {"chebi_id": CHEBI_ATP, "name": "ATP"},
                {
                    "chebi_id": CHEBI_COA,
                    "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS",
                    "name": "CoA",
                },
            ],
            right_participants=[
                {
                    "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                    "name": "acetyl-CoA",
                },
                {"chebi_id": CHEBI_AMP, "name": "AMP"},
                {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_detects_adp_phosphate_coupling(self):
        classifier = AcidThiolLigase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "O=C([O-])CCC(=O)[O-]", "name": "succinate"},
                {
                    "chebi_id": CHEBI_COA,
                    "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS",
                    "name": "CoA",
                },
                {
                    "chebi_id": CHEBI_ATP,
                    "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                    "name": "ATP",
                },
            ],
            right_participants=[
                {
                    "smiles": "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCC(=O)[O-]",
                    "name": "succinyl-CoA",
                },
                {
                    "chebi_id": "CHEBI:456216",
                    "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
                    "name": "ADP",
                },
                {"chebi_id": CHEBI_PHOSPHATE, "smiles": "O=P([O-])([O-])O", "name": "phosphate"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is True

    def test_rejects_without_thioester_product(self):
        classifier = AcidThiolLigase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC(=O)[O-]", "name": "acetate"},
                {"chebi_id": CHEBI_ATP, "name": "ATP"},
                {"chebi_id": CHEBI_COA, "name": "CoA"},
            ],
            right_participants=[
                {"smiles": "CC(=O)[O-]", "name": "acetate"},
                {"chebi_id": CHEBI_AMP, "name": "AMP"},
                {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_thioester_hydrolysis(self):
        classifier = CarboxylicEsterHydrolase()
        reaction = make_reaction(
            left_participants=[
                {"smiles": "CC(=O)SCCN", "name": "thioester"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {"smiles": "CC(=O)[O-]", "name": "acetate"},
                {"smiles": "NCCS", "name": "aminoethanethiol"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_bisphosphate_monoester_hydrolysis(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "OCC1OC(O)(COP(=O)([O-])[O-])C(O)C(O)C1OP(=O)([O-])[O-]",
                    "name": "fructose 1,6-bisphosphate",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "OCC1OC(O)C(O)C(O)C1OP(=O)([O-])[O-]",
                    "name": "fructose 6-phosphate",
                },
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_mixed_hydrolysis_with_ammonia_release(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "NC(=O)N1C=C[C@H](CO[P](=O)([O-])O[P](=O)([O-])O[P](=O)([O-])[O-])O1",
                    "name": "dCTP-like substrate",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "O=C1NCC(COP(=O)([O-])[O-])O1",
                    "name": "dephosphorylated nucleotide",
                },
                {
                    "chebi_id": CHEBI_DIPHOSPHATE,
                    "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                    "name": "diphosphate",
                },
                {"smiles": "[NH4+]", "name": "ammonium"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_nadp_phosphatase_like_dephosphorylation(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](OP(=O)([O-])[O-])[C@@H]3O)[C@@H](O)[C@H]2O)=C1",
                    "name": "NADP",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)=C1",
                    "name": "NAD",
                },
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_multiple_disconnected_phosphoryl_domains(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])O)[C@@H](OP(=O)([O-])OP(=O)([O-])[O-])[C@H]2O)C(=O)N1",
                    "name": "guanosine 3',5'-bis(diphosphate)",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1",
                    "name": "GDP",
                },
                {
                    "chebi_id": CHEBI_DIPHOSPHATE,
                    "smiles": "O=P([O-])([O-])OP(=O)([O-])O",
                    "name": "diphosphate",
                },
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False

    def test_rejects_multiwater_hydrolysis_beyond_anhydride_cleavage(self):
        classifier = Pyrophosphatase()
        reaction = make_reaction(
            left_participants=[
                {
                    "smiles": "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1",
                    "name": "GTP",
                },
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
                {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            ],
            right_participants=[
                {
                    "smiles": "[H]C(=O)NC1=C(N[C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)N=C(N)NC1=O",
                    "name": "formylamino phosphoribosyl pyrimidinone",
                },
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
                {
                    "chebi_id": CHEBI_PHOSPHATE,
                    "smiles": "O=P([O-])([O-])O",
                    "name": "phosphate",
                },
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
                {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
            ],
        )

        result = classifier.check_membership_impl(reaction)

        assert result.is_member is False
