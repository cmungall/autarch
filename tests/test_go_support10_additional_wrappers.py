"""Focused positive tests for additional GO-backed wrappers with support >= 10."""

from __future__ import annotations

from typing import Callable

import pytest
from typing_extensions import TypedDict

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H2O2,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.abc_type_transporter_activity import ABCTypeTransporterActivity
from autarch.ontology.active_monoatomic_ion_transmembrane_transporter_activity import (
    ActiveMonoatomicIonTransmembraneTransporterActivity,
)
from autarch.ontology.active_transmembrane_transporter_activity import (
    ActiveTransmembraneTransporterActivity,
)
from autarch.ontology.alcohol_dehydrogenase_nad_or_nadp import (
    AlcoholDehydrogenaseNADOrNADP,
)
from autarch.ontology.aldehyde_dehydrogenase_nad_or_nadp import (
    AldehydeDehydrogenaseNADOrNADP,
)
from autarch.ontology.antioxidant_activity import AntioxidantActivity
from autarch.ontology.atp_dependent_activity import ATPDependentActivity
from autarch.ontology.carbohydrate_kinase_activity import CarbohydrateKinaseActivity
from autarch.ontology.carbohydrate_transmembrane_transporter_activity import (
    CarbohydrateTransmembraneTransporterActivity,
)
from autarch.ontology.catalytic_activity_acting_on_a_rrna import (
    CatalyticActivityActingOnARRNA,
)
from autarch.ontology.catalytic_activity_acting_on_a_nucleic_acid import (
    CatalyticActivityActingOnANucleicAcid,
)
from autarch.ontology.catalytic_activity_acting_on_a_protein import (
    CatalyticActivityActingOnAProtein,
)
from autarch.ontology.catalytic_activity_acting_on_a_trna import (
    CatalyticActivityActingOnATRNA,
)
from autarch.ontology.catalytic_activity_acting_on_rna import (
    CatalyticActivityActingOnRNA,
)
from autarch.ontology.deacylase_activity import DeacylaseActivity
from autarch.ontology.ligase_activity_forming_carbon_oxygen_bonds import (
    LigaseActivityFormingCarbonOxygenBonds,
)
from autarch.ontology.monoatomic_cation_transmembrane_transporter_activity import (
    MonoatomicCationTransmembraneTransporterActivity,
)
from autarch.ontology.monoatomic_ion_transmembrane_transporter_activity import (
    MonoatomicIonTransmembraneTransporterActivity,
)
from autarch.ontology.n_acetyltransferase import NAcetyltransferase
from autarch.ontology.n_acyltransferase import NAcyltransferase
from autarch.ontology.nucleoside_triphosphate_diphosphatase_activity import (
    NucleosideTriphosphateDiphosphataseActivity,
)
from autarch.ontology.o_acetyltransferase import OAcetyltransferase
from autarch.ontology.o_acyltransferase import OAcyltransferase
from autarch.ontology.proton_transmembrane_transporter_activity import (
    ProtonTransmembraneTransporterActivity,
)
from autarch.ontology.rna_dihydrouridine_synthase import RNADihydrouridineSynthase
from autarch.ontology.rrna_methyltransferase import RRNAMethyltransferase
from autarch.ontology.s_adenosylmethionine_dependent_methyltransferase import (
    SAdenosylmethionineDependentMethyltransferase,
)
from autarch.ontology.s_methyltransferase import SMethyltransferase
from autarch.ontology.lipase import Lipase
from autarch.ontology.transmembrane_transporter_activity import (
    TransmembraneTransporterActivity,
)
from autarch.ontology.transporter_activity import TransporterActivity
from autarch.ontology.trna_dihydrouridine_synthase import TRNADihydrouridineSynthase
from autarch.ontology.trna_methyltransferase import TRNAMethyltransferase

CHEBI_POLYPHOSPHATE = "CHEBI:16838"
CHEBI_D_XYLOSE = "CHEBI:53455"
CHEBI_PEPTIDYL_PROLINE_CIS = "CHEBI:83155"
CHEBI_PEPTIDYL_PROLINE_TRANS = "CHEBI:83154"
CHEBI_CYTIDINE_IN_TRNA = "CHEBI:82748"
CHEBI_METHYLCYTIDINE_IN_TRNA = "CHEBI:74483"
CHEBI_TRNA_DIHYDROURIDINE = "CHEBI:74443"
CHEBI_TRNA_URIDINE = "CHEBI:65315"
CHEBI_D_TYROSYL_TRNA = "CHEBI:78723"
CHEBI_D_TYROSINE = "CHEBI:58570"
CHEBI_TRNA = "CHEBI:78442"
CHEBI_D_GLUCOSE = "CHEBI:4167"
CHEBI_D_GLUCOSE_6_PHOSPHATE = "CHEBI:61548"
CHEBI_ACETYL_COA = "CHEBI:57288"
CHEBI_NAT_A_SUBSTRATE = "CHEBI:64739"
CHEBI_NAT_A_PRODUCT = "CHEBI:133375"
CHEBI_METHIONINE = "CHEBI:57844"
CHEBI_S_METHYL_METHIONINE = "CHEBI:58252"
CHEBI_PSEUDOURIDINE_IN_RRNA = "CHEBI:65314"
CHEBI_METHYLPSEUDOURIDINE_IN_RRNA = "CHEBI:74890"
CHEBI_SUPEROXIDE = "CHEBI:18421"
CHEBI_TYROSINE = "CHEBI:58315"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    name: str
    count: int
    location: str
    polymer_index: str
    polymer_type: PolymerType


def _reaction(left: list[ParticipantInput], right: list[ParticipantInput], label: str | None = None) -> Reaction:
    kwargs = {}
    if label is not None:
        kwargs["label"] = label
    return Reaction(
        left_participants=[Participant(**participant) for participant in left],
        right_participants=[Participant(**participant) for participant in right],
        **kwargs,
    )


def _atp_coupled_transport(substrate_left: ParticipantInput, substrate_right: ParticipantInput) -> Reaction:
    return _reaction(
        left=[
            substrate_left,
            {"chebi_id": CHEBI_ATP, "smiles": "P", "name": "ATP"},
            {"chebi_id": CHEBI_H2O, "smiles": "[H]O[H]", "name": "water"},
        ],
        right=[
            substrate_right,
            {"chebi_id": CHEBI_ADP, "smiles": "P", "name": "ADP"},
            {"chebi_id": CHEBI_POLYPHOSPHATE, "name": "polyphosphate"},
            {"chebi_id": CHEBI_H_PLUS, "smiles": "[H+]", "name": "hydron"},
        ],
    )


def _protein_isomerase_reaction() -> Reaction:
    return _reaction(
        left=[{"chebi_id": CHEBI_PEPTIDYL_PROLINE_CIS}],
        right=[{"chebi_id": CHEBI_PEPTIDYL_PROLINE_TRANS}],
    )


def _trna_methylation_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_CYTIDINE_IN_TRNA},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_METHYLCYTIDINE_IN_TRNA},
            {"chebi_id": CHEBI_SAH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _abc_carbohydrate_transport_reaction() -> Reaction:
    return _atp_coupled_transport(
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "out",
        },
        {
            "chebi_id": CHEBI_D_XYLOSE,
            "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)C1",
            "name": "D-xylopyranose",
            "location": "in",
        },
    )


def _hydron_translocation_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": "CHEBI:29034"},
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": "CHEBI:29033"},
            {"chebi_id": CHEBI_H2O},
        ],
    )


def _deacylase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_D_TYROSYL_TRNA, "name": "D-tyrosyl-tRNA", "polymer_type": PolymerType.TRNA},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"chebi_id": CHEBI_D_TYROSINE, "name": "D-tyrosine"},
            {"chebi_id": CHEBI_TRNA, "name": "tRNA", "polymer_type": PolymerType.TRNA},
        ],
    )


def _alcohol_dehydrogenase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"smiles": "CCO", "name": "ethanol"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )


def _antioxidant_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_SUPEROXIDE},
            {"chebi_id": CHEBI_H_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_O2},
            {"chebi_id": CHEBI_H2O2},
        ],
    )


def _aminoacyl_trna_ligase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_TRNA, "polymer_type": PolymerType.TRNA, "name": "tRNA"},
            {"chebi_id": CHEBI_TYROSINE, "name": "L-tyrosine"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": CHEBI_TRNA, "polymer_type": PolymerType.TRNA, "name": "tyrosyl-tRNA"},
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
        label="tRNA + L-tyrosine + ATP = tyrosyl-tRNA + AMP + diphosphate",
    )


def _s_methyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_METHIONINE},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_S_METHYL_METHIONINE},
            {"chebi_id": CHEBI_SAH},
        ],
    )


def _polysialic_o_acetyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "name": "sialic acid polymer",
                "polymer_type": PolymerType.SIALIC_ACID_POLYMER,
                "polymer_index": "n",
            },
            {"chebi_id": CHEBI_ACETYL_COA, "name": "acetyl-CoA"},
        ],
        right=[
            {
                "name": "acetylated sialic acid polymer",
                "polymer_type": PolymerType.SIALIC_ACID_POLYMER,
                "polymer_index": "n",
            },
            {"chebi_id": CHEBI_COA, "name": "CoA"},
        ],
    )


def _rrna_methylation_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_PSEUDOURIDINE_IN_RRNA},
            {"chebi_id": CHEBI_SAM},
        ],
        right=[
            {"chebi_id": CHEBI_METHYLPSEUDOURIDINE_IN_RRNA},
            {"chebi_id": CHEBI_SAH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _lipase_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": "CHEBI:17855",
                "smiles": "CC(=O)OCC(COC(C)=O)OC(C)=O",
                "name": "triacetin",
            },
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
            {"chebi_id": CHEBI_H2O, "smiles": "O", "name": "water"},
        ],
        right=[
            {"chebi_id": "CHEBI:15734", "smiles": "OCC(O)CO", "name": "glycerol"},
            {"chebi_id": "CHEBI:18059", "smiles": "CC(=O)O", "name": "acetate", "count": 3},
        ],
    )


def _rna_dihydrouridine_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_TRNA_DIHYDROURIDINE},
            {"chebi_id": CHEBI_NADP_PLUS},
        ],
        right=[
            {"chebi_id": CHEBI_TRNA_URIDINE},
            {"chebi_id": CHEBI_NADPH},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _nucleoside_triphosphate_diphosphatase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_ATP},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[
            {"chebi_id": CHEBI_AMP},
            {"chebi_id": CHEBI_DIPHOSPHATE},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _n_acetyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_NAT_A_SUBSTRATE, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_ACETYL_COA},
        ],
        right=[
            {"chebi_id": CHEBI_NAT_A_PRODUCT, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_COA},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _carbohydrate_kinase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_D_GLUCOSE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_D_GLUCOSE_6_PHOSPHATE},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _aldehyde_dehydrogenase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"smiles": "CC=O", "name": "acetaldehyde"},
            {"smiles": "OP(=O)([O-])[O-]", "name": "phosphate"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
        right=[
            {"smiles": "CC(=O)OP(=O)([O-])[O-]", "name": "acetyl phosphate"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )


@pytest.mark.parametrize(
    ("cls", "reaction_factory"),
    [
        (CatalyticActivityActingOnAProtein, _protein_isomerase_reaction),
        (SAdenosylmethionineDependentMethyltransferase, _trna_methylation_reaction),
        (CatalyticActivityActingOnANucleicAcid, _trna_methylation_reaction),
        (CatalyticActivityActingOnRNA, _trna_methylation_reaction),
        (CatalyticActivityActingOnARRNA, _rrna_methylation_reaction),
        (CatalyticActivityActingOnATRNA, _trna_methylation_reaction),
        (TransporterActivity, _abc_carbohydrate_transport_reaction),
        (TransmembraneTransporterActivity, _abc_carbohydrate_transport_reaction),
        (ActiveTransmembraneTransporterActivity, _abc_carbohydrate_transport_reaction),
        (ATPDependentActivity, _abc_carbohydrate_transport_reaction),
        (MonoatomicIonTransmembraneTransporterActivity, _hydron_translocation_reaction),
        (MonoatomicCationTransmembraneTransporterActivity, _hydron_translocation_reaction),
        (DeacylaseActivity, _deacylase_reaction),
        (ABCTypeTransporterActivity, _abc_carbohydrate_transport_reaction),
        (ProtonTransmembraneTransporterActivity, _hydron_translocation_reaction),
        (TRNAMethyltransferase, _trna_methylation_reaction),
        (AlcoholDehydrogenaseNADOrNADP, _alcohol_dehydrogenase_reaction),
        (AntioxidantActivity, _antioxidant_reaction),
        (LigaseActivityFormingCarbonOxygenBonds, _aminoacyl_trna_ligase_reaction),
        (ActiveMonoatomicIonTransmembraneTransporterActivity, _hydron_translocation_reaction),
        (CarbohydrateTransmembraneTransporterActivity, _abc_carbohydrate_transport_reaction),
        (SMethyltransferase, _s_methyltransferase_reaction),
        (OAcyltransferase, _polysialic_o_acetyltransferase_reaction),
        (OAcetyltransferase, _polysialic_o_acetyltransferase_reaction),
        (RRNAMethyltransferase, _rrna_methylation_reaction),
        (RNADihydrouridineSynthase, _rna_dihydrouridine_reaction),
        (TRNADihydrouridineSynthase, _rna_dihydrouridine_reaction),
        (NucleosideTriphosphateDiphosphataseActivity, _nucleoside_triphosphate_diphosphatase_reaction),
        (NAcyltransferase, _n_acetyltransferase_reaction),
        (NAcetyltransferase, _n_acetyltransferase_reaction),
        (CarbohydrateKinaseActivity, _carbohydrate_kinase_reaction),
        (AldehydeDehydrogenaseNADOrNADP, _aldehyde_dehydrogenase_reaction),
        (Lipase, _lipase_reaction),
    ],
)
def test_go_support10_additional_wrapper_positive_examples(
    cls: type,
    reaction_factory: Callable[[], Reaction],
) -> None:
    result = cls().check_membership(reaction_factory())
    assert result.is_member is True
