"""Focused tests for the next five grouping classifiers."""

from abc import abstractmethod
from typing import Any

from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_ATP,
    CHEBI_CMP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_NH4,
)
from autarch.ontology.hydrolase_acting_on_carbon_nitrogen_but_not_peptide_bonds_in_cyclic_amidines import (
    HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmidines,
)
from autarch.ontology.hydrolase_acting_on_carbon_nitrogen_but_not_peptide_bonds_in_nitriles import (
    HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInNitriles,
)
from autarch.ontology.nucleotidyltransferase import Nucleotidyltransferase
from autarch.ontology.oxidoreductase_acting_on_nadh_or_nadph_quinone_or_similar_compound_as_acceptor import (
    OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor,
)
from autarch.ontology.transferring_other_glycosyl_groups import (
    TransferringOtherGlycosylGroups,
)
from autarch.ontology.reaction import ReactionClass
from autarch.validation.go_term_validator import GoTermValidator

ATP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])"
    "[C@@H](O)[C@H]1O"
)
ADP_SMILES = (
    "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
)
GTP_SMILES = (
    "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)"
    "C(=O)N1"
)
GDP_SMILES = (
    "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C(=O)N1"
)
NADH_SMILES = (
    "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)"
    "[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1"
)
NAD_PLUS_SMILES = (
    "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)"
    "[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)=C1"
)
NADPH_SMILES = (
    "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)"
    "[C@H](O)[C@@H]3OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=CC1"
)
NADP_PLUS_SMILES = (
    "NC(=O)C1=CC=C[N+]([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)"
    "[C@H](O)[C@@H]3OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)=C1"
)
DIPHOSPHATE_SMILES = "O=P([O-])([O-])OP(=O)([O-])O"
PHOSPHATE_SMILES = "O=P([O-])([O-])[O-]"
CMP_SMILES = "NC1=NC(=O)N([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]2O)C=C1"


class ParticipantInput(TypedDict, total=False):
    """Subset of participant fields used in classifier tests."""

    chebi_id: str
    smiles: str
    name: str
    count: int
    polymer_index: str
    polymer_type: PolymerType
    monomer: str


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput]) -> Reaction:
    """Build a reaction from lightweight participant dictionaries."""
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
    )


def test_ec_165_positive_quinone_reduction() -> None:
    cls = OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=C1C=CC(=O)C([O-])=C1",
                "name": "2-hydroxy-1,4-benzoquinone",
            },
            {"chebi_id": CHEBI_NADH, "smiles": NADH_SMILES, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"smiles": "OC1=CC(O)=C(O)C=C1", "name": "benzene-1,2,4-triol"},
            {"chebi_id": CHEBI_NAD_PLUS, "smiles": NAD_PLUS_SMILES, "name": "NAD(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "quinone" in result.explanation.lower() or "similar acceptor" in result.explanation.lower()


def test_ec_165_positive_similar_acceptor_reduction() -> None:
    cls = OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:59513",
                "smiles": "[H][C@]1([C@@H](O)CO)OC(=O)C([O-])=C1[O]",
                "name": "monodehydro-L-ascorbate",
            },
            {"chebi_id": CHEBI_NADH, "smiles": NADH_SMILES, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "chebi_id": "CHEBI:38290",
                "smiles": "[H][C@]1([C@@H](O)CO)OC(=O)C(O)=C1[O-]",
                "name": "L-ascorbate",
            },
            {"chebi_id": CHEBI_NAD_PLUS, "smiles": NAD_PLUS_SMILES, "name": "NAD(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_165_rejects_generic_aldehyde_reduction() -> None:
    cls = OxidoreductaseActingOnNADHOrNADPHQuinoneOrSimilarCompoundAsAcceptor()
    reaction = _reaction(
        left=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_NADPH, "smiles": NADPH_SMILES, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_NADP_PLUS, "smiles": NADP_PLUS_SMILES, "name": "NADP(+)"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_354_positive_cyclic_amidine_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmidines()
    reaction = _reaction(
        left=[
            {"smiles": "CN1CC(=O)NC1=N", "name": "creatinine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {"smiles": "CN1CC(=O)NC1=O", "name": "N-methylhydantoin"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "cyclic amidine" in result.explanation.lower()


def test_ec_354_positive_nucleoside_deamination() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmidines()
    reaction = _reaction(
        left=[
            {
                "smiles": "NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
                "name": "adenosine",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "smiles": "O=C1NC=NC2=C1N=CN2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
                "name": "inosine",
            },
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_354_rejects_linear_amide_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInCyclicAmidines()
    reaction = _reaction(
        left=[
            {"smiles": "CC(N)=O", "name": "acetamide"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "CC(=O)[O-]", "name": "acetate"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_355_positive_nitrile_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInNitriles()
    reaction = _reaction(
        left=[
            {"smiles": "N#CC[C@H]([NH3+])C(=O)[O-]", "name": "3-cyano-L-alanine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water", "count": 2},
        ],
        right=[
            {"smiles": "[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]", "name": "L-aspartate"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "nitrile" in result.explanation.lower()


def test_ec_355_positive_thiocyanate_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInNitriles()
    reaction = _reaction(
        left=[
            {"smiles": "[S-]C#N", "name": "thiocyanate"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron", "count": 2},
        ],
        right=[
            {"smiles": "O=C=S", "name": "carbonyl sulfide"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_355_rejects_amide_hydrolysis() -> None:
    cls = HydrolaseActingOnCarbonNitrogenButNotPeptideBondsInNitriles()
    reaction = _reaction(
        left=[
            {"smiles": "NC(=O)CC(=O)[O-]", "name": "asparagine"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"smiles": "[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]", "name": "aspartate"},
            {"chebi_id": CHEBI_NH4, "smiles": "[H][N+]([H])([H])[H]", "name": "ammonium"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_277_positive_nucleotidyl_transfer_to_sugar_phosphate() -> None:
    cls = Nucleotidyltransferase()
    reaction = _reaction(
        left=[
            {
                "smiles": "O=P([O-])([O-])O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                "name": "alpha-D-glucose 1-phosphate",
            },
            {"chebi_id": CHEBI_GTP, "smiles": GTP_SMILES, "name": "GTP"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
        right=[
            {
                "smiles": "NC1=NC2=C(/N=C\\N2[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OC3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)"
                "[C@@H](O)[C@H]2O)C(=O)N1",
                "name": "GDP-D-glucose",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": DIPHOSPHATE_SMILES, "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "nucleotidyl" in result.explanation.lower()


def test_ec_277_positive_polymerase_like_branch() -> None:
    cls = Nucleotidyltransferase()
    reaction = _reaction(
        left=[
            {
                "name": "RNA",
                "polymer_type": PolymerType.RNA,
                "polymer_index": "n",
                "monomer": "ribonucleotide",
            },
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
        ],
        right=[
            {
                "name": "RNA",
                "polymer_type": PolymerType.RNA,
                "polymer_index": "n+1",
                "monomer": "ribonucleotide",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "smiles": DIPHOSPHATE_SMILES, "name": "diphosphate"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_277_rejects_simple_kinase() -> None:
    cls = Nucleotidyltransferase()
    reaction = _reaction(
        left=[
            {"smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "name": "glucose"},
            {"chebi_id": CHEBI_ATP, "smiles": ATP_SMILES, "name": "ATP"},
        ],
        right=[
            {
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                "name": "glucose 6-phosphate",
            },
            {"chebi_id": CHEBI_ADP, "smiles": ADP_SMILES, "name": "ADP"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


def test_ec_2499_positive_kdo_transfer() -> None:
    cls = TransferringOtherGlycosylGroups()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:58603",
                "smiles": "CCCCCCCCCCC[C@@H](O)CC(=O)N[C@H]1[C@H](OC[C@H]2O[C@H](OP(=O)([O-])[O-])[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)"
                "[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO)[C@@H](OP(=O)([O-])[O-])[C@@H]1OC(=O)C[C@H](O)CCCCCCCCCCC",
                "name": "lipid IVA",
            },
            {
                "chebi_id": "CHEBI:85987",
                "smiles": "[H][C@]1([C@H](O)CO)O[C@](OP(=O)([O-])OC[C@H]2O[C@@H](N3C=CC(N)=NC3=O)[C@H](O)[C@@H]2O)(C(=O)[O-])C[C@@H](O)[C@H]1O",
                "name": "CMP-Kdo",
            },
        ],
        right=[
            {
                "chebi_id": "CHEBI:60364",
                "smiles": "CCCCCCCCCCC[C@@H](O)CC(=O)N[C@H]1[C@H](OC[C@H]2O[C@H](OP(=O)([O-])[O-])[C@H](NC(=O)C[C@H](O)CCCCCCCCCCC)"
                "[C@@H](OC(=O)C[C@H](O)CCCCCCCCCCC)[C@@H]2O)O[C@H](CO[C@]2(C(=O)[O-])C[C@@H](O)[C@@H](O)[C@@H]([C@H](O)CO)O2)"
                "[C@@H](OP(=O)([O-])[O-])[C@@H]1OC(=O)C[C@H](O)CCCCCCCCCCC",
                "name": "Kdo-lipid IVA",
            },
            {"chebi_id": CHEBI_CMP, "smiles": CMP_SMILES, "name": "CMP"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True
    assert "glycosyl" in result.explanation.lower()


def test_ec_2499_positive_lipid_linked_oligosaccharide_transfer() -> None:
    cls = TransferringOtherGlycosylGroups()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:57570",
                "name": "dolichyl diphosphooligosaccharide",
            },
            {
                "chebi_id": "CHEBI:50347",
                "name": "L-asparaginyl-[protein]",
                "polymer_type": PolymerType.PROTEIN,
                "monomer": "amino acid",
            },
        ],
        right=[
            {
                "chebi_id": "CHEBI:132529",
                "name": "glycosyl-asparaginyl-[protein]",
                "polymer_type": PolymerType.PROTEIN,
                "monomer": "amino acid",
            },
            {
                "chebi_id": "CHEBI:57497",
                "name": "dolichyl diphosphate",
            },
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is True


def test_ec_2499_rejects_hexosyltransferase_chemistry() -> None:
    cls = TransferringOtherGlycosylGroups()
    reaction = _reaction(
        left=[
            {
                "chebi_id": "CHEBI:18066",
                "smiles": "O=P([O-])([O-])OC[C@H]1O[C@@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O",
                "name": "UDP-D-glucose",
            },
            {"smiles": "CO", "name": "methanol"},
        ],
        right=[
            {"chebi_id": "CHEBI:58223", "name": "UDP"},
            {"smiles": "CO[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O", "name": "methyl glucoside"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )

    result = cls.check_membership_impl(reaction)

    assert result.is_member is False


class FakeOntologyAdapter:
    """Simple adapter that can return labels for GO or EC terms."""

    def __init__(self, label_map: dict[str, str]):
        self._label_map = label_map

    def label(self, term_id: str) -> str | None:
        return self._label_map.get(term_id)


def make_ec_only_reaction_class(
    class_name: str,
    module_name: str,
    ec_prefix: str,
    docstring: str,
) -> type[ReactionClass]:
    """Create a minimal concrete EC-only class for validator tests."""

    def check_membership_impl(self: ReactionClass, reaction: Reaction) -> Any:
        raise AssertionError("not called")

    attrs: dict[str, Any] = {
        "__module__": module_name,
        "__doc__": docstring,
        "GO_ID": None,
        "EC_NUMBER_PREFIX": ec_prefix,
        "check_membership_impl": abstractmethod(check_membership_impl),
    }
    return type(class_name, (ReactionClass,), attrs)


def test_validator_accepts_ec_only_alignment_when_go_id_is_absent() -> None:
    cls = make_ec_only_reaction_class(
        class_name="TransferringOtherGlycosylGroups",
        module_name="synthetic.transferring_other_glycosyl_groups",
        ec_prefix="2.4.99.-",
        docstring="Transferring other glycosyl groups.",
    )
    validator = GoTermValidator(
        adapter=FakeOntologyAdapter({"EC:2.4.99": "Transferring other glycosyl groups"}),
    )

    result = validator.validate_class(cls)

    assert result.matches
    assert result.file_matches_class
    assert result.description_matches_class
