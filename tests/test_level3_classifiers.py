"""
Tests for EC level 3 classifiers: CarboxylicEsterHydrolase, OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor, Aminopeptidase.
"""

import pytest

from autarch.ontology.carboxylic_ester_hydrolase import CarboxylicEsterHydrolase
from autarch.ontology.oxidoreductase_acting_on_the_ch_oh_group_of_donors_oxygen_as_acceptor import (
    OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor,
)
from autarch.ontology.aminopeptidase import Aminopeptidase
from autarch.datamodel import Reaction, Participant


def create_reaction(left_participants, right_participants):
    """Helper to create reaction objects."""
    return Reaction(
        left_participants=[
            Participant(chebi_id=p.get("chebi_id"), name=p.get("name"), smiles=p.get("smiles"))
            for p in left_participants
        ],
        right_participants=[
            Participant(chebi_id=p.get("chebi_id"), name=p.get("name"), smiles=p.get("smiles"))
            for p in right_participants
        ]
    )


class TestCarboxylicEsterHydrolase:
    """Test CarboxylicEsterHydrolase classifier."""

    def test_simple_ester_hydrolysis(self):
        """Test aromatic ester hydrolysis with explicit structures."""
        esterase = CarboxylicEsterHydrolase()
        reaction = create_reaction(
            left_participants=[
                {"chebi_id": "CHEBI:15377", "name": "H2O", "smiles": "O"},
                {"name": "methyl salicylate", "smiles": "COC(=O)c1ccccc1O"},
            ],
            right_participants=[
                {"name": "methanol", "smiles": "CO"},
                {"name": "salicylate", "smiles": "O=C([O-])c1ccccc1O"},
                {"chebi_id": "CHEBI:15378", "name": "H(+)", "smiles": "[H+]"},
            ],
        )

        result = esterase.check_membership_impl(reaction)
        assert result.is_member is True
        assert "ester" in result.explanation.lower()

    def test_lactone_hydrolysis(self):
        """Test lactone ring opening."""
        esterase = CarboxylicEsterHydrolase()
        reaction = create_reaction(
            left_participants=[
                {"chebi_id": "CHEBI:15377", "name": "H2O", "smiles": "O"},
                {"name": "3,4-dihydrocoumarin", "smiles": "O=C1OCCCc2ccccc21"},
            ],
            right_participants=[
                {"name": "3-(2-hydroxyphenyl)propanoate", "smiles": "O=C([O-])CCCc1ccccc1O"},
                {"chebi_id": "CHEBI:15378", "name": "H(+)", "smiles": "[H+]"},
            ],
        )

        result = esterase.check_membership_impl(reaction)
        assert result.is_member is True

    @pytest.mark.skip(reason="Name-based detection removed - classifiers use ChEBI IDs only")
    def test_acetate_ester_by_name(self):
        """Test ester detection by name pattern."""
        esterase = CarboxylicEsterHydrolase()
        reaction = create_reaction(
            left_participants=[
                {"name": "H2O"},
                {"name": "ethyl acetate"}
            ],
            right_participants=[
                {"name": "acetate"},
                {"name": "ethanol"}
            ]
        )

        result = esterase.check_membership_impl(reaction)
        assert result.is_member is True

    def test_non_ester_reaction(self):
        """Test that non-ester reactions are rejected.""" 
        esterase = CarboxylicEsterHydrolase()
        reaction = create_reaction(
            left_participants=[
                {"name": "glucose", "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
                {"chebi_id": "CHEBI:15422", "name": "ATP"}
            ],
            right_participants=[
                {"name": "glucose-6-phosphate", "smiles": "O=P([O-])([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
                {"chebi_id": "CHEBI:16761", "name": "ADP"},
            ],
        )

        result = esterase.check_membership_impl(reaction)
        assert result.is_member is False

    def test_no_water_reactant(self):
        """Test rejection when no water present."""
        esterase = CarboxylicEsterHydrolase()
        reaction = create_reaction(
            left_participants=[
                {"name": "ethyl acetate", "smiles": "CCOC(C)=O"}
            ],
            right_participants=[
                {"name": "acetate", "smiles": "CC(=O)[O-]"},
                {"name": "ethanol", "smiles": "CCO"},
            ],
        )

        result = esterase.check_membership_impl(reaction)
        assert result.is_member is False
        assert "no water" in result.explanation.lower()


class TestOxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor:
    """Test OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor classifier."""

    def test_glucose_oxidase(self):
        """Test RHEA:10036 - glucose oxidation."""
        alcohol_oxidase = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = create_reaction(
            left_participants=[
                {"chebi_id": "CHEBI:15379", "name": "O2"},
                {"chebi_id": "CHEBI:17925", "name": "D-glucose"}
            ],
            right_participants=[
                {"chebi_id": "CHEBI:16240", "name": "H2O2"},
                {"chebi_id": "CHEBI:16217", "name": "D-glucono-1,5-lactone"}
            ]
        )
        
        result = alcohol_oxidase.check_membership_impl(reaction)
        assert result.is_member is True
        assert "alcohol oxidase" in result.explanation.lower()

    def test_vanillyl_alcohol_oxidation(self):
        """Test vanillyl alcohol → vanillin."""
        alcohol_oxidase = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = create_reaction(
            left_participants=[
                {"chebi_id": "CHEBI:15379", "name": "O2"},
                {"chebi_id": "CHEBI:18353", "name": "4-hydroxy-3-methoxy-benzenemethanol"}
            ],
            right_participants=[
                {"chebi_id": "CHEBI:16240", "name": "H2O2"},
                {"chebi_id": "CHEBI:18346", "name": "vanillin"}
            ]
        )
        
        result = alcohol_oxidase.check_membership_impl(reaction)
        assert result.is_member is True

    @pytest.mark.skip(reason="Name-based detection removed - classifiers use ChEBI IDs only")
    def test_methanol_oxidation_by_name(self):
        """Test methanol oxidation by name patterns."""
        alcohol_oxidase = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = create_reaction(
            left_participants=[
                {"name": "O2"},
                {"name": "methanol"}
            ],
            right_participants=[
                {"name": "H2O2"},
                {"name": "formaldehyde"}
            ]
        )

        result = alcohol_oxidase.check_membership_impl(reaction)
        assert result.is_member is True

    def test_no_oxygen_reactant(self):
        """Test rejection when no O2 present."""
        alcohol_oxidase = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = create_reaction(
            left_participants=[
                {"name": "glucose"},
                {"name": "NAD+"}
            ],
            right_participants=[
                {"name": "gluconate"},
                {"name": "NADH"}
            ]
        )
        
        result = alcohol_oxidase.check_membership_impl(reaction)
        assert result.is_member is False
        assert "no oxygen" in result.explanation.lower()

    def test_no_peroxide_product(self):
        """Test rejection when no H2O2 produced (using ChEBI IDs)."""
        alcohol_oxidase = OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor()
        reaction = create_reaction(
            left_participants=[
                {"chebi_id": "CHEBI:15379", "name": "O2"},  # O2
                {"chebi_id": "CHEBI:17234", "name": "glucose"}
            ],
            right_participants=[
                {"chebi_id": "CHEBI:15377", "name": "water"},  # H2O, not H2O2
                {"chebi_id": "CHEBI:16658", "name": "gluconate"}
            ]
        )

        result = alcohol_oxidase.check_membership_impl(reaction)
        assert result.is_member is False
        assert "hydrogen peroxide" in result.explanation.lower()


class TestAminopeptidase:
    """Test Aminopeptidase classifier."""

    @pytest.mark.skip(reason="Polymer detection requires name-based patterns - classifiers use ChEBI IDs only")
    def test_tryptophanyl_aminopeptidase(self):
        """Test RHEA:72999 - N-terminal tryptophan cleavage."""
        aminopeptidase = Aminopeptidase()
        reaction = create_reaction(
            left_participants=[
                {"name": "H2O"},
                {"name": "an N-terminal L-tryptophanyl-L-alpha-aminoacyl-[peptide]"}
            ],
            right_participants=[
                {"name": "an N-terminal L-alpha-aminoacyl-[peptide]"},
                {"chebi_id": "CHEBI:16828", "name": "L-tryptophan"}
            ]
        )

        result = aminopeptidase.check_membership_impl(reaction)
        assert result.is_member is True
        assert "aminopeptidase" in result.explanation.lower()

    @pytest.mark.skip(reason="Polymer detection requires name-based patterns - classifiers use ChEBI IDs only")
    def test_arginyl_aminopeptidase(self):
        """Test RHEA:78707 - N-terminal arginine cleavage."""
        aminopeptidase = Aminopeptidase()
        reaction = create_reaction(
            left_participants=[
                {"name": "H2O"},
                {"name": "an N-terminal L-arginyl-L-aminoacyl-[protein]"}
            ],
            right_participants=[
                {"name": "an N-terminal L-alpha-aminoacyl-[protein]"},
                {"chebi_id": "CHEBI:16467", "name": "L-arginine"}
            ]
        )

        result = aminopeptidase.check_membership_impl(reaction)
        assert result.is_member is True

    @pytest.mark.skip(reason="Polymer detection requires name-based patterns - classifiers use ChEBI IDs only")
    def test_peptide_amino_acid_cleavage(self):
        """Test general peptide cleavage pattern."""
        aminopeptidase = Aminopeptidase()
        reaction = create_reaction(
            left_participants=[
                {"name": "H2O"},
                {"name": "dipeptide"}
            ],
            right_participants=[
                {"name": "peptide"},
                {"name": "alanine"}
            ]
        )

        result = aminopeptidase.check_membership_impl(reaction)
        assert result.is_member is True

    def test_non_peptidase_reaction(self):
        """Test rejection of non-peptidase reactions."""
        aminopeptidase = Aminopeptidase()
        reaction = create_reaction(
            left_participants=[
                {"name": "H2O"},
                {"name": "glucose"}
            ],
            right_participants=[
                {"name": "fructose"},
                {"name": "galactose"}
            ]
        )
        
        result = aminopeptidase.check_membership_impl(reaction)
        assert result.is_member is False

    def test_no_water_reactant(self):
        """Test rejection when no water present."""
        aminopeptidase = Aminopeptidase()
        reaction = create_reaction(
            left_participants=[
                {"name": "peptide"}
            ],
            right_participants=[
                {"name": "amino acid"},
                {"name": "shorter peptide"}
            ]
        )
        
        result = aminopeptidase.check_membership_impl(reaction)
        assert result.is_member is False
        assert "no water" in result.explanation.lower()
