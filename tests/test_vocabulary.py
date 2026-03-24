"""Tests for the chemical vocabulary package."""


class TestVocabularyImports:
    """Test that vocabulary can be imported correctly."""

    def test_import_from_main_vocabulary(self):
        """Test importing common molecules from main vocabulary."""
        from autarch.vocabulary import ATP, ADP, NAD, NADH, glucose

        assert ATP.chebi_id == "CHEBI:15422"
        assert ATP.name == "ATP"
        assert ADP.chebi_id == "CHEBI:16761"
        assert NAD.chebi_id == "CHEBI:57540"
        assert NADH.chebi_id == "CHEBI:57945"
        assert glucose.chebi_id == "CHEBI:17234"

    def test_import_from_energy_module(self):
        """Test importing from energy submodule."""
        from autarch.vocabulary.energy import ATP, GTP, CTP, UTP, Pi, PPi

        assert ATP.chebi_id == "CHEBI:15422"
        assert GTP.chebi_id == "CHEBI:15996"
        assert CTP.chebi_id == "CHEBI:17677"
        assert UTP.chebi_id == "CHEBI:15713"
        assert Pi.chebi_id == "CHEBI:43474"
        assert PPi.chebi_id == "CHEBI:33019"

    def test_import_from_cofactors_module(self):
        """Test importing from cofactors submodule."""
        from autarch.vocabulary.cofactors import (
            NAD,
            NADH,
            FAD,
            FADH2,
            CoA,
            acetyl_CoA,
            SAM,
            SAH,
        )

        assert NAD.chebi_id == "CHEBI:57540"
        assert NADH.chebi_id == "CHEBI:57945"
        assert FAD.chebi_id == "CHEBI:57692"
        assert FADH2.chebi_id == "CHEBI:57618"
        assert CoA.chebi_id == "CHEBI:15346"
        assert acetyl_CoA.chebi_id == "CHEBI:15351"
        assert SAM.chebi_id == "CHEBI:15414"
        assert SAH.chebi_id == "CHEBI:16680"

    def test_import_from_common_module(self):
        """Test importing from common submodule."""
        from autarch.vocabulary.common import H2O, H, OH, O2, CO2, NH3

        assert H2O.chebi_id == "CHEBI:15377"
        assert H.chebi_id == "CHEBI:15378"
        assert OH.chebi_id == "CHEBI:16234"
        assert O2.chebi_id == "CHEBI:15379"
        assert CO2.chebi_id == "CHEBI:16526"
        assert NH3.chebi_id == "CHEBI:16134"

    def test_import_from_metabolites_module(self):
        """Test importing from metabolites submodule."""
        from autarch.vocabulary.metabolites import (
            glucose,
            pyruvate,
            lactate,
            citrate,
            succinate,
            fumarate,
            glutamate,
            glutamine,
        )

        assert glucose.chebi_id == "CHEBI:17234"
        assert pyruvate.chebi_id == "CHEBI:15361"
        assert lactate.chebi_id == "CHEBI:24996"
        assert citrate.chebi_id == "CHEBI:16947"
        assert succinate.chebi_id == "CHEBI:15741"
        assert fumarate.chebi_id == "CHEBI:18012"
        assert glutamate.chebi_id == "CHEBI:16015"
        assert glutamine.chebi_id == "CHEBI:18050"

    def test_import_from_ions_module(self):
        """Test importing from ions submodule."""
        from autarch.vocabulary.ions import (
            Na,
            K,
            Ca,
            Mg,
            Fe2,
            Fe3,
            Zn2,
            Cl,
            phosphate,
            sulfate,
        )

        assert Na.chebi_id == "CHEBI:29101"
        assert K.chebi_id == "CHEBI:29103"
        assert Ca.chebi_id == "CHEBI:29108"
        assert Mg.chebi_id == "CHEBI:18420"
        assert Fe2.chebi_id == "CHEBI:29033"
        assert Fe3.chebi_id == "CHEBI:29034"
        assert Zn2.chebi_id == "CHEBI:29105"
        assert Cl.chebi_id == "CHEBI:17996"
        assert phosphate.chebi_id == "CHEBI:43474"
        assert sulfate.chebi_id == "CHEBI:16189"

    def test_all_molecules_dictionary(self):
        """Test the ALL_MOLECULES dictionary."""
        from autarch.vocabulary import ALL_MOLECULES

        # Check it's a dictionary
        assert isinstance(ALL_MOLECULES, dict)

        # Check some key molecules are present
        assert "ATP" in ALL_MOLECULES
        assert "NAD" in ALL_MOLECULES
        assert "glucose" in ALL_MOLECULES
        assert "Na" in ALL_MOLECULES

        # Check values are Molecule objects
        atp = ALL_MOLECULES["ATP"]
        assert atp.chebi_id == "CHEBI:15422"
        assert atp.name == "ATP"

    def test_backward_compatibility(self):
        """Test that molecules are still available from formula module."""
        from autarch.formula import ATP, ADP, NAD, NADH, glucose

        # These should still work for backward compatibility
        assert ATP.chebi_id == "CHEBI:15422"
        assert ADP.chebi_id == "CHEBI:16761"
        assert NAD.chebi_id == "CHEBI:57540"
        assert NADH.chebi_id == "CHEBI:57945"
        assert glucose.chebi_id == "CHEBI:17234"


class TestVocabularyUsage:
    """Test using vocabulary molecules in reactions."""

    def test_reaction_with_vocabulary(self):
        """Test creating reactions with vocabulary molecules."""
        from autarch.vocabulary import ATP, ADP, Pi, H2O

        # ATP hydrolysis
        reaction = (ATP + H2O >> ADP + Pi).build()

        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 2
        assert reaction.left_participants[0].chebi_id == "CHEBI:15422"
        assert reaction.left_participants[1].chebi_id == "CHEBI:15377"
        assert reaction.right_participants[0].chebi_id == "CHEBI:16761"
        assert reaction.right_participants[1].chebi_id == "CHEBI:43474"

    def test_complex_pathway_with_vocabulary(self):
        """Test metabolic pathway with vocabulary."""
        from autarch.vocabulary.metabolites import glucose, glucose_6P, fructose_6P
        from autarch.vocabulary import ATP, ADP, H

        # Glycolysis steps
        hexokinase = glucose + ATP >> glucose_6P + ADP + H
        isomerase = glucose_6P >> fructose_6P

        r1 = hexokinase.build()
        r2 = isomerase.build()

        assert r1.left_participants[0].name == "glucose"
        assert r1.right_participants[0].name == "G6P"

        assert r2.left_participants[0].name == "G6P"
        assert r2.right_participants[0].name == "F6P"

    def test_transport_with_ions(self):
        """Test transport reactions with ions."""
        from autarch.vocabulary.ions import Na, K
        from autarch.vocabulary import ATP, ADP, Pi

        # Na/K pump
        pump = (Na["in"] * 3 + K["out"] * 2 + ATP) >> (
            Na["out"] * 3 + K["in"] * 2 + ADP + Pi
        )

        reaction = pump.build()

        # Check stoichiometry
        na_in = [p for p in reaction.left_participants if p.chebi_id == "CHEBI:29101"]
        assert len(na_in) == 1
        assert na_in[0].count == 3
        assert na_in[0].location == "in"

        k_out = [p for p in reaction.left_participants if p.chebi_id == "CHEBI:29103"]
        assert len(k_out) == 1
        assert k_out[0].count == 2
        assert k_out[0].location == "out"
