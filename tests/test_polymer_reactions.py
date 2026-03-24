"""Tests for polymer reaction parsing and handling.

Tests the (n) stoichiometry notation for polymer reactions like:
- RNA(n) + ribonucleoside 5'-diphosphate => RNA(n+1) + phosphate
- [(1->4)-alpha-D-glucosyl](n) + H2O => [(1->4)-alpha-D-glucosyl](n-1) + D-glucose
"""

import pytest
from autarch.etl.label_parser import (
    parse_reaction_label,
    parse_stoichiometry_from_term,
    parse_polymer_index,
    infer_polymer_info,
    has_polymer_notation,
)
from autarch.etl.rhea_tsv_etl import RheaTSVETL
from autarch.datamodel import PolymerType


class TestPolymerNotationParsing:
    """Test parsing of (n), (n+1), (n-1) notation."""

    def test_has_polymer_notation(self):
        """Test detection of polymer notation in labels."""
        assert has_polymer_notation("RNA(n) + ATP => RNA(n+1) + ADP")
        assert has_polymer_notation("glucan(n) + H2O => glucan(n-1) + glucose")
        assert has_polymer_notation("n ATP + H2O => n ADP + n phosphate")
        assert not has_polymer_notation("ATP + H2O => ADP + phosphate")

    def test_parse_polymer_index_simple(self):
        """Test parsing simple polymer index from molecule names.

        parse_polymer_index returns (base_name, polymer_index) tuple.
        """
        base, idx = parse_polymer_index("RNA(n)")
        assert base == "RNA"
        assert idx == "n"

        base, idx = parse_polymer_index("RNA(n+1)")
        assert base == "RNA"
        assert idx == "n+1"

        base, idx = parse_polymer_index("RNA(n-1)")
        assert base == "RNA"
        assert idx == "n-1"

        base, idx = parse_polymer_index("polysulfur(n-2)")
        assert base == "polysulfur"
        assert idx == "n-2"

        base, idx = parse_polymer_index("ATP")
        assert base == "ATP"
        assert idx is None

    def test_parse_stoichiometry_coefficient(self):
        """Test parsing variable stoichiometry like 'n ATP' or '2n H2O'."""
        coef, name = parse_stoichiometry_from_term("n ATP")
        assert coef == "n"
        assert name == "ATP"

        coef, name = parse_stoichiometry_from_term("2n H2O")
        assert coef == "2n"
        assert name == "H2O"

        coef, name = parse_stoichiometry_from_term("(n-1) sulfur")
        assert coef == "n-1"
        assert name == "sulfur"

        coef, name = parse_stoichiometry_from_term("ATP")
        assert coef is None
        assert name == "ATP"


class TestPolymerTypeInference:
    """Test inference of polymer type from molecule names."""

    def test_infer_rna(self):
        """Test RNA polymer type inference."""
        ptype, linkage, monomer = infer_polymer_info("RNA")
        assert ptype == "rna"
        assert monomer == "ribonucleotide"

    def test_infer_trna(self):
        """Test tRNA polymer type inference."""
        ptype, linkage, monomer = infer_polymer_info("tRNA")
        assert ptype == "trna"
        assert monomer == "ribonucleotide"

    def test_infer_glucan_with_linkage(self):
        """Test glucan with linkage notation."""
        ptype, linkage, monomer = infer_polymer_info("[(1->4)-alpha-D-glucosyl]")
        assert ptype == "glucan"
        assert linkage == "1->4"
        assert "glucosyl" in monomer.lower()

    def test_infer_chitin(self):
        """Test chitin polymer type inference."""
        ptype, linkage, monomer = infer_polymer_info("[(1->4)-N-acetyl-beta-D-glucosaminyl]")
        assert ptype == "chitin"
        assert linkage == "1->4"

    def test_infer_sialic_acid(self):
        """Test sialic acid polymer inference."""
        ptype, linkage, monomer = infer_polymer_info("[N-acetyl-alpha-D-neuraminosyl-(2->8)]")
        assert ptype == "sialic_acid_polymer"
        assert linkage == "2->8"


class TestReactionLabelParsing:
    """Test parsing of full reaction labels.

    parse_reaction_label returns (left_participants, right_participants) tuple.
    """

    def test_parse_rna_polymerization(self):
        """Test parsing RNA polymerization reaction - RHEA:22098."""
        label = "RNA(n) + a ribonucleoside 5'-diphosphate => RNA(n+1) + phosphate"
        left, right = parse_reaction_label(label)

        assert len(left) == 2
        assert len(right) == 2

        # Check RNA(n) on left
        rna_left = next(p for p in left if "RNA" in p.name)
        assert rna_left.polymer_index == "n"

        # Check RNA(n+1) on right
        rna_right = next(p for p in right if "RNA" in p.name)
        assert rna_right.polymer_index == "n+1"

    def test_parse_glucan_hydrolysis(self):
        """Test parsing glucan hydrolysis reaction."""
        label = "[(1->4)-alpha-D-glucosyl](n) + H2O = [(1->4)-alpha-D-glucosyl](n-1) + alpha-D-glucose"
        left, right = parse_reaction_label(label)

        assert len(left) == 2
        assert len(right) == 2

        # Check polymer indices
        glucan_left = next(p for p in left if "glucosyl" in p.name.lower())
        assert glucan_left.polymer_index == "n"

        glucan_right = next(p for p in right if "glucosyl" in p.name.lower())
        assert glucan_right.polymer_index == "n-1"

    def test_parse_variable_stoichiometry(self):
        """Test parsing reactions with n ATP, 2n H2O etc."""
        label = "n ATP + protein => n ADP + phospho-protein"
        left, right = parse_reaction_label(label)

        atp = next(p for p in left if "ATP" in p.name)
        assert atp.stoichiometry == "n"

        adp = next(p for p in right if "ADP" in p.name)
        assert adp.stoichiometry == "n"


class TestPolymerReactionLoading:
    """Test loading polymer reactions from RHEA data."""

    @pytest.fixture
    def etl(self, tmp_path):
        """Create ETL instance with test cache."""
        return RheaTSVETL(cache_dir=tmp_path / "rhea_tsv")

    def test_load_polymer_labels(self, etl):
        """Test that polymer labels are loaded correctly."""
        labels = etl.get_polymer_labels()
        # Should find polymer reactions
        assert len(labels) > 0

        # Check for known polymer reaction
        # RHEA:22098 is RNA(n) + ribonucleoside => RNA(n+1) + phosphate
        if "22098" in labels or "22097" in labels:
            # Found the RNA polymerase reaction
            label = labels.get("22098") or labels.get("22097")
            assert "RNA" in label
            assert "(n)" in label or "(n+1)" in label

    def test_reaction_from_label(self, etl):
        """Test creating Reaction objects from labels."""
        label = "RNA(n) + a ribonucleoside 5'-diphosphate = RNA(n+1) + phosphate"
        reaction = etl.reaction_from_label("22098", label)

        assert reaction is not None
        assert len(reaction.left_participants) == 2
        assert len(reaction.right_participants) == 2

        # Check polymer info is populated
        rna_left = next(
            (p for p in reaction.left_participants if p.name and "RNA" in p.name),
            None
        )
        assert rna_left is not None
        assert rna_left.polymer_index == "n"
        assert rna_left.polymer_type == PolymerType.RNA


class TestSpecificPolymerReactions:
    """Test specific polymer reactions from RHEA."""

    def test_rhea_22098_rna_polymerization(self):
        """Test RHEA:22098 - polynucleotide phosphorylase reaction.

        RNA(n) + ribonucleoside 5'-diphosphate => RNA(n+1) + phosphate
        EC: 2.7.7.8
        """
        from autarch.etl.label_parser import parse_reaction_label

        label = "RNA(n) + a ribonucleoside 5'-diphosphate = RNA(n+1) + phosphate"
        left, right = parse_reaction_label(label)

        # Left side: RNA(n) and ribonucleoside
        assert len(left) == 2
        rna_n = next(p for p in left if "RNA" in p.name)
        assert rna_n.polymer_index == "n"

        # Right side: RNA(n+1) and phosphate
        assert len(right) == 2
        rna_n1 = next(p for p in right if "RNA" in p.name)
        assert rna_n1.polymer_index == "n+1"

        # This shows polymer growth: n -> n+1
        assert rna_n.polymer_index == "n"
        assert rna_n1.polymer_index == "n+1"

    def test_rhea_10257_glucan_kinase(self):
        """Test RHEA:10257 - glucan kinase reaction.

        [(1->4)-6-phospho-alpha-D-glucosyl](n) + n ATP + n H2O =>
        [(1->4)-3,6-bisphospho-alpha-D-glucosyl](n) + n AMP + n phosphate + 2n H(+)
        EC: 2.7.9.5
        """
        from autarch.etl.label_parser import parse_reaction_label

        label = (
            "[(1->4)-6-phospho-alpha-D-glucosyl](n) + n ATP + n H2O = "
            "[(1->4)-3,6-bisphospho-alpha-D-glucosyl](n) + n AMP + n phosphate + 2n H(+)"
        )
        left, right = parse_reaction_label(label)

        # Check left side has polymer and variable stoichiometry
        glucan_left = next(p for p in left if "glucosyl" in p.name.lower())
        assert glucan_left.polymer_index == "n"

        atp = next(p for p in left if "ATP" in p.name)
        assert atp.stoichiometry == "n"

        # Check right side
        glucan_right = next(p for p in right if "glucosyl" in p.name.lower())
        assert glucan_right.polymer_index == "n"

        h_plus = next(p for p in right if "H(+)" in p.name)
        assert h_plus.stoichiometry == "2n"
