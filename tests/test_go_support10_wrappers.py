"""Focused tests for GO-backed wrappers with transitive support above 10."""

from __future__ import annotations

from typing import Callable
from typing_extensions import TypedDict

import pytest

from autarch.datamodel import Participant, PolymerType, Reaction
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_COA,
    CHEBI_DIPHOSPHATE,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NADH,
    CHEBI_NADPH,
    CHEBI_NADP_PLUS,
    CHEBI_NAD_PLUS,
    CHEBI_O2,
    CHEBI_SAH,
    CHEBI_SAM,
)
from autarch.ontology.acetyltransferase import Acetyltransferase
from autarch.ontology.adenylyltransferase import Adenylyltransferase
from autarch.ontology.carbon_nitrogen_lyase import CarbonNitrogenLyase
from autarch.ontology.carbon_oxygen_lyase_acting_on_phosphates import (
    CarbonOxygenLyaseActingOnPhosphates,
)
from autarch.ontology.coa_ligase import CoALigase
from autarch.ontology.dioxygenase import Dioxygenase
from autarch.ontology.hydrolase_acting_on_carbon_carbon_bonds import (
    HydrolaseActingOnCarbonCarbonBonds,
)
from autarch.ontology.hydrolase_acting_on_ether_bonds import HydrolaseActingOnEtherBonds
from autarch.ontology.intramolecular_oxidoreductase_transposing_c_c_bonds import (
    IntramolecularOxidoreductaseTransposingCCBonds,
)
from autarch.ontology.intramolecular_transferase import IntramolecularTransferase
from autarch.ontology.n_methyltransferase import NMethyltransferase
from autarch.ontology.nucleobase_containing_compound_kinase import (
    NucleobaseContainingCompoundKinase,
)
from autarch.ontology.o_methyltransferase import OMethyltransferase
from autarch.ontology.oxidoreductase_acting_on_nadh_or_nadph import (
    OxidoreductaseActingOnNADHOrNADPH,
)
from autarch.ontology.phosphoric_ester_hydrolase import PhosphoricEsterHydrolase
from autarch.ontology.rna_methyltransferase import RNAMethyltransferase

CHEBI_NITROBENZENE = "CHEBI:27798"
CHEBI_CATECHOL = "CHEBI:18135"
CHEBI_NITRITE = "CHEBI:16301"
CHEBI_GPP = "CHEBI:58057"
CHEBI_ALPHA_PINENE_MINUS = "CHEBI:28660"
CHEBI_EPOXIDE = "CHEBI:35762"
CHEBI_DIOL = "CHEBI:23824"
CHEBI_FUMARYLACETOACETATE = "CHEBI:18034"
CHEBI_PHENYLALANINE = "CHEBI:17295"
CHEBI_NH3 = "CHEBI:16134"
CHEBI_PHOSPHOETHANOLAMINE = "CHEBI:58190"
CHEBI_N_METHYL_PHOSPHOETHANOLAMINE = "CHEBI:57781"
CHEBI_THYMIDINE = "CHEBI:17748"
CHEBI_DTMP = "CHEBI:63528"
CHEBI_MENAQUINONE = "CHEBI:16374"
CHEBI_MENAQUINOL = "CHEBI:18151"
CHEBI_CROTONYL_COA = "CHEBI:26348"
CHEBI_CYTIDINE_IN_TRNA = "CHEBI:82748"
CHEBI_METHYLCYTIDINE_IN_TRNA = "CHEBI:74483"
CHEBI_NAT_A_SUBSTRATE = "CHEBI:64739"
CHEBI_NAT_A_PRODUCT = "CHEBI:133375"
CHEBI_CHOLATE = "CHEBI:29747"
CHEBI_CHOLOYL_COA = "CHEBI:57373"
CHEBI_GAMMA_TOCOPHEROL = "CHEBI:18185"
CHEBI_ALPHA_TOCOPHEROL = "CHEBI:18145"


class ParticipantInput(TypedDict, total=False):
    chebi_id: str
    smiles: str
    name: str
    count: int
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


def _dioxygenase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_NITROBENZENE, "name": "nitrobenzene"},
            {"chebi_id": CHEBI_NADH, "name": "NADH"},
            {"chebi_id": CHEBI_O2, "name": "dioxygen"},
        ],
        right=[
            {"chebi_id": CHEBI_CATECHOL, "name": "catechol"},
            {"chebi_id": CHEBI_NITRITE, "name": "nitrite"},
            {"chebi_id": CHEBI_NAD_PLUS, "name": "NAD+"},
        ],
    )


def _phosphoric_ester_hydrolase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"smiles": "C=CCOP(=O)(O)OP(=O)(O)O", "name": "allyl diphosphate"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"smiles": "C=CCO", "name": "allyl alcohol"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


def _carbon_oxygen_lyase_acting_on_phosphates_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_GPP,
                "smiles": "CC(C)=CCC/C(C)=C/COP(=O)([O-])OP(=O)([O-])[O-]",
                "name": "geranyl diphosphate",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_ALPHA_PINENE_MINUS,
                "smiles": "CC1=CC[C@H]2C[C@@H]1C2(C)C",
                "name": "alpha-pinene",
            },
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


def _intramolecular_transferase_reaction() -> Reaction:
    return _reaction(
        left=[{"name": "glucose 1-phosphate"}],
        right=[{"name": "glucose 6-phosphate"}],
        label="glucose 1-phosphate = glucose 6-phosphate",
    )


def _n_methyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_PHOSPHOETHANOLAMINE, "name": "phosphoethanolamine"},
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {"chebi_id": CHEBI_N_METHYL_PHOSPHOETHANOLAMINE, "name": "N-methylethanolamine phosphate"},
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )


def _carbon_nitrogen_lyase_reaction() -> Reaction:
    return _reaction(
        left=[{"chebi_id": CHEBI_PHENYLALANINE, "name": "L-phenylalanine"}],
        right=[{"name": "trans-cinnamate"}, {"chebi_id": CHEBI_NH3, "name": "ammonia"}],
        label="L-phenylalanine = trans-cinnamate + ammonia",
    )


def _nucleobase_containing_compound_kinase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_THYMIDINE},
            {"chebi_id": CHEBI_ATP},
        ],
        right=[
            {"chebi_id": CHEBI_DTMP},
            {"chebi_id": CHEBI_ADP},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _oxidoreductase_acting_on_nadh_or_nadph_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_MENAQUINONE, "smiles": "CC1=C(C)C(=O)C=CC1=O", "name": "menaquinone"},
            {"chebi_id": CHEBI_NADPH, "name": "NADPH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
        right=[
            {"chebi_id": CHEBI_MENAQUINOL, "smiles": "CC1=C(C)C(O)=CC=C1O", "name": "menaquinol"},
            {"chebi_id": CHEBI_NADP_PLUS, "name": "NADP+"},
        ],
    )


def _hydrolase_acting_on_carbon_carbon_bonds_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_FUMARYLACETOACETATE},
            {"chebi_id": CHEBI_H2O},
        ],
        right=[{"name": "fumarate"}, {"name": "acetoacetate"}],
        label="fumarylacetoacetate + water = fumarate + acetoacetate",
    )


def _intramolecular_oxidoreductase_transposing_c_c_bonds_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "smiles": "CCCCCCCC/C=C\\CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "(3Z)-dodecenoyl-CoA",
            },
        ],
        right=[
            {
                "chebi_id": CHEBI_CROTONYL_COA,
                "smiles": "CCCCCCCCC/C=C/C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](N2C=NC3=C2N=CN=C3N)[C@H](O)[C@@H]1OP(=O)([O-])[O-]",
                "name": "(2E)-dodecenoyl-CoA",
            },
        ],
    )


def _hydrolase_acting_on_ether_bonds_reaction() -> Reaction:
    return _reaction(
        left=[
            {"smiles": "C1OC1c1ccccc1", "name": "styrene oxide"},
            {"chebi_id": CHEBI_H2O, "name": "water"},
        ],
        right=[
            {"smiles": "OC(CO)c1ccccc1", "name": "phenylethane-1,2-diol"},
        ],
    )


def _adenylyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": "CHEBI:46858", "name": "L-tyrosyl-[protein]"},
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
        ],
        right=[
            {"chebi_id": "CHEBI:83624", "name": "O-(5'-adenylyl)-L-tyrosyl-[protein]"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


def _acetyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {"chebi_id": CHEBI_NAT_A_SUBSTRATE, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": "CHEBI:57288"},
        ],
        right=[
            {"chebi_id": CHEBI_NAT_A_PRODUCT, "polymer_type": PolymerType.PROTEIN},
            {"chebi_id": CHEBI_COA},
            {"chebi_id": CHEBI_H_PLUS},
        ],
    )


def _o_methyltransferase_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_GAMMA_TOCOPHEROL,
                "smiles": "CC1=C(O)\\C=C2\\CC[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O\\C2=C\\1C",
                "name": "gamma-tocopherol",
            },
            {"chebi_id": CHEBI_SAM, "name": "SAM"},
        ],
        right=[
            {
                "chebi_id": CHEBI_ALPHA_TOCOPHEROL,
                "smiles": "CC1=C2\\CC[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O\\C2=C(C)\\C(C)=C\\1O",
                "name": "alpha-tocopherol",
            },
            {"chebi_id": CHEBI_SAH, "name": "SAH"},
            {"chebi_id": CHEBI_H_PLUS, "name": "hydron"},
        ],
    )


def _rna_methyltransferase_reaction() -> Reaction:
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


def _coa_ligase_reaction() -> Reaction:
    return _reaction(
        left=[
            {
                "chebi_id": CHEBI_CHOLATE,
                "smiles": "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@@]3(C)[C@@]([H])(CC[C@]3([H])[C@H](C)CCC(=O)[O-])[C@]1([H])[C@H](O)C2",
                "name": "cholate",
            },
            {"chebi_id": CHEBI_ATP, "name": "ATP"},
            {"chebi_id": CHEBI_COA, "name": "CoA"},
        ],
        right=[
            {
                "chebi_id": CHEBI_CHOLOYL_COA,
                "smiles": "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@@]3(C)[C@@]([H])(CC[C@]3([H])[C@]([H])(C)CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](N4C=NC5=C4N=CN=C5N)[C@H](O)[C@@H]3OP(=O)([O-])[O-])[C@]1([H])[C@H](O)C2",
                "name": "choloyl-CoA",
            },
            {"chebi_id": CHEBI_AMP, "name": "AMP"},
            {"chebi_id": CHEBI_DIPHOSPHATE, "name": "diphosphate"},
        ],
    )


@pytest.mark.parametrize(
    ("cls", "reaction_factory"),
    [
        (Dioxygenase, _dioxygenase_reaction),
        (PhosphoricEsterHydrolase, _phosphoric_ester_hydrolase_reaction),
        (CarbonOxygenLyaseActingOnPhosphates, _carbon_oxygen_lyase_acting_on_phosphates_reaction),
        (IntramolecularTransferase, _intramolecular_transferase_reaction),
        (NMethyltransferase, _n_methyltransferase_reaction),
        (CarbonNitrogenLyase, _carbon_nitrogen_lyase_reaction),
        (NucleobaseContainingCompoundKinase, _nucleobase_containing_compound_kinase_reaction),
        (OxidoreductaseActingOnNADHOrNADPH, _oxidoreductase_acting_on_nadh_or_nadph_reaction),
        (HydrolaseActingOnCarbonCarbonBonds, _hydrolase_acting_on_carbon_carbon_bonds_reaction),
        (IntramolecularOxidoreductaseTransposingCCBonds, _intramolecular_oxidoreductase_transposing_c_c_bonds_reaction),
        (HydrolaseActingOnEtherBonds, _hydrolase_acting_on_ether_bonds_reaction),
        (Adenylyltransferase, _adenylyltransferase_reaction),
        (Acetyltransferase, _acetyltransferase_reaction),
        (OMethyltransferase, _o_methyltransferase_reaction),
        (RNAMethyltransferase, _rna_methyltransferase_reaction),
        (CoALigase, _coa_ligase_reaction),
    ],
)
def test_go_support10_wrapper_positive_examples(
    cls: type,
    reaction_factory: Callable[[], Reaction],
) -> None:
    reaction = reaction_factory()
    result = cls().check_membership(reaction)
    assert result.is_member is True
