"""Tests for ChEBI normalization module."""

import gzip
import tempfile
from pathlib import Path


from autarch.etl.chebi_normalization import (
    COMMON_MOLECULE_CHEBI,
    EquationChEBILookup,
    normalize_chebi_by_name,
    parse_equation_chebi_ids,
    parse_equation_to_position_map,
)


class TestNormalizeChebiByName:
    """Tests for name-based ChEBI lookup."""

    def test_proton_variants(self):
        """Test various proton/hydron names."""
        assert normalize_chebi_by_name("H(+)") == "CHEBI:15378"
        assert normalize_chebi_by_name("H+") == "CHEBI:15378"
        assert normalize_chebi_by_name("proton") == "CHEBI:15378"
        assert normalize_chebi_by_name("hydron") == "CHEBI:15378"
        assert normalize_chebi_by_name("hydrogen ion") == "CHEBI:15378"

    def test_water_variants(self):
        """Test water name variants."""
        assert normalize_chebi_by_name("H2O") == "CHEBI:15377"
        assert normalize_chebi_by_name("water") == "CHEBI:15377"
        assert normalize_chebi_by_name("h(2)o") == "CHEBI:15377"

    def test_phosphate_variants(self):
        """Test phosphate species."""
        assert normalize_chebi_by_name("phosphate") == "CHEBI:43474"
        assert normalize_chebi_by_name("Pi") == "CHEBI:43474"
        assert normalize_chebi_by_name("diphosphate") == "CHEBI:33019"
        assert normalize_chebi_by_name("PPi") == "CHEBI:33019"

    def test_nucleotides_atp_family(self):
        """Test ATP/ADP/AMP lookup."""
        assert normalize_chebi_by_name("ATP") == "CHEBI:30616"
        assert normalize_chebi_by_name("ADP") == "CHEBI:456216"
        assert normalize_chebi_by_name("AMP") == "CHEBI:456215"
        assert normalize_chebi_by_name("adenosine triphosphate") == "CHEBI:30616"

    def test_nucleotides_gtp_family(self):
        """Test GTP/GDP/GMP lookup."""
        assert normalize_chebi_by_name("GTP") == "CHEBI:37565"
        assert normalize_chebi_by_name("GDP") == "CHEBI:58189"
        assert normalize_chebi_by_name("GMP") == "CHEBI:58115"

    def test_nucleotides_utp_family(self):
        """Test UTP/UDP/UMP lookup."""
        assert normalize_chebi_by_name("UTP") == "CHEBI:46398"
        assert normalize_chebi_by_name("UDP") == "CHEBI:58223"
        assert normalize_chebi_by_name("UMP") == "CHEBI:57865"

    def test_nad_cofactors(self):
        """Test NAD/NADH/NADP/NADPH lookup."""
        assert normalize_chebi_by_name("NAD(+)") == "CHEBI:57540"
        assert normalize_chebi_by_name("NAD+") == "CHEBI:57540"
        assert normalize_chebi_by_name("NADH") == "CHEBI:57945"
        assert normalize_chebi_by_name("NADP(+)") == "CHEBI:58349"
        assert normalize_chebi_by_name("NADPH") == "CHEBI:57783"

    def test_fad_cofactors(self):
        """Test FAD/FMN cofactor lookup."""
        assert normalize_chebi_by_name("FAD") == "CHEBI:57692"
        assert normalize_chebi_by_name("FADH2") == "CHEBI:58307"
        assert normalize_chebi_by_name("FMN") == "CHEBI:58210"
        assert normalize_chebi_by_name("FMNH2") == "CHEBI:57618"

    def test_coa_variants(self):
        """Test CoA and derivatives."""
        assert normalize_chebi_by_name("CoA") == "CHEBI:57287"
        assert normalize_chebi_by_name("coenzyme A") == "CHEBI:57287"
        assert normalize_chebi_by_name("acetyl-CoA") == "CHEBI:57288"
        assert normalize_chebi_by_name("malonyl-CoA") == "CHEBI:57384"

    def test_gases(self):
        """Test common gases."""
        assert normalize_chebi_by_name("O2") == "CHEBI:15379"
        assert normalize_chebi_by_name("oxygen") == "CHEBI:15379"
        assert normalize_chebi_by_name("CO2") == "CHEBI:16526"
        assert normalize_chebi_by_name("H2") == "CHEBI:18276"
        assert normalize_chebi_by_name("NH3") == "CHEBI:16134"
        assert normalize_chebi_by_name("NH4(+)") == "CHEBI:28938"

    def test_sulfur_compounds(self):
        """Test sulfur species."""
        assert normalize_chebi_by_name("sulfur") == "CHEBI:17909"
        assert normalize_chebi_by_name("S") == "CHEBI:17909"
        assert normalize_chebi_by_name("hydrogen sulfide") == "CHEBI:16136"
        assert normalize_chebi_by_name("sulfate") == "CHEBI:16189"

    def test_common_metabolites(self):
        """Test common metabolic intermediates."""
        assert normalize_chebi_by_name("acetate") == "CHEBI:30089"
        assert normalize_chebi_by_name("pyruvate") == "CHEBI:15361"
        assert normalize_chebi_by_name("citrate") == "CHEBI:16947"
        assert normalize_chebi_by_name("2-oxoglutarate") == "CHEBI:16810"

    def test_case_insensitive(self):
        """Test that lookup is case-insensitive."""
        assert normalize_chebi_by_name("atp") == "CHEBI:30616"
        assert normalize_chebi_by_name("ATP") == "CHEBI:30616"
        assert normalize_chebi_by_name("Atp") == "CHEBI:30616"
        assert normalize_chebi_by_name("WATER") == "CHEBI:15377"

    def test_whitespace_handling(self):
        """Test whitespace is stripped."""
        assert normalize_chebi_by_name("  ATP  ") == "CHEBI:30616"
        assert normalize_chebi_by_name("\twater\n") == "CHEBI:15377"

    def test_unknown_molecule(self):
        """Test unknown molecules return None."""
        assert normalize_chebi_by_name("unknown_molecule_xyz") is None
        assert normalize_chebi_by_name("this_does_not_exist") is None

    def test_empty_input(self):
        """Test empty/None input handling."""
        assert normalize_chebi_by_name("") is None
        assert normalize_chebi_by_name(None) is None
        assert normalize_chebi_by_name("   ") is None


class TestParseEquationChebiIds:
    """Tests for EQUATION parsing."""

    def test_simple_equation(self):
        """Test simple A = B equation."""
        left, right = parse_equation_chebi_ids("CHEBI:16459 = CHEBI:31011")
        assert left == ["CHEBI:16459"]
        assert right == ["CHEBI:31011"]

    def test_multiple_participants(self):
        """Test equation with multiple participants on each side."""
        left, right = parse_equation_chebi_ids(
            "CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938"
        )
        assert left == ["CHEBI:16459", "CHEBI:15377"]
        assert right == ["CHEBI:31011", "CHEBI:28938"]

    def test_arrow_directions(self):
        """Test different arrow types."""
        # Undefined (=)
        left, right = parse_equation_chebi_ids("CHEBI:111 = CHEBI:222")
        assert left == ["CHEBI:111"]
        assert right == ["CHEBI:222"]

        # Left-to-right (=>)
        left, right = parse_equation_chebi_ids("CHEBI:111 => CHEBI:222")
        assert left == ["CHEBI:111"]
        assert right == ["CHEBI:222"]

        # Bidirectional (<=>)
        left, right = parse_equation_chebi_ids("CHEBI:111 <=> CHEBI:222")
        assert left == ["CHEBI:111"]
        assert right == ["CHEBI:222"]

    def test_stoichiometry_prefix(self):
        """Test numeric stoichiometry prefix."""
        left, right = parse_equation_chebi_ids("2 CHEBI:15378 = CHEBI:31011")
        assert left == ["CHEBI:15378", "CHEBI:15378"]
        assert right == ["CHEBI:31011"]

        left, right = parse_equation_chebi_ids(
            "CHEBI:111 + 3 CHEBI:222 = 2 CHEBI:333"
        )
        assert left == ["CHEBI:111", "CHEBI:222", "CHEBI:222", "CHEBI:222"]
        assert right == ["CHEBI:333", "CHEBI:333"]

    def test_variable_stoichiometry(self):
        """Test variable stoichiometry (n, 2n, etc.)."""
        left, right = parse_equation_chebi_ids("n CHEBI:15378 = CHEBI:31011")
        assert left == ["CHEBI:15378"]  # Variable stoich -> single entry
        assert right == ["CHEBI:31011"]

        left, right = parse_equation_chebi_ids("2n CHEBI:15378 = CHEBI:31011")
        assert left == ["CHEBI:15378"]
        assert right == ["CHEBI:31011"]

    def test_comma_separated(self):
        """Test comma-separated ChEBI IDs."""
        left, right = parse_equation_chebi_ids(
            "CHEBI:29950,CHEBI:29950 = CHEBI:50058"
        )
        assert left == ["CHEBI:29950", "CHEBI:29950"]
        assert right == ["CHEBI:50058"]

    def test_complex_equation(self):
        """Test complex real-world equation."""
        equation = (
            "CHEBI:57384 + CHEBI:57288 + 2 CHEBI:57783 + 2 CHEBI:15378 "
            "<=> CHEBI:57560 + CHEBI:16526 + 2 CHEBI:58349 + 6 CHEBI:57287"
        )
        left, right = parse_equation_chebi_ids(equation)
        assert len(left) == 6  # 1 + 1 + 2 + 2
        assert "CHEBI:57384" in left
        assert "CHEBI:57288" in left
        assert left.count("CHEBI:57783") == 2
        assert left.count("CHEBI:15378") == 2

        assert len(right) == 10  # 1 + 1 + 2 + 6
        assert right.count("CHEBI:57287") == 6

    def test_empty_input(self):
        """Test empty input handling."""
        assert parse_equation_chebi_ids("") == ([], [])
        assert parse_equation_chebi_ids("   ") == ([], [])

    def test_invalid_format(self):
        """Test invalid equation format."""
        assert parse_equation_chebi_ids("no equals sign") == ([], [])
        assert parse_equation_chebi_ids("just text") == ([], [])

    def test_whitespace_handling(self):
        """Test various whitespace patterns."""
        left, right = parse_equation_chebi_ids(
            "  CHEBI:111  +  CHEBI:222  =  CHEBI:333  "
        )
        assert left == ["CHEBI:111", "CHEBI:222"]
        assert right == ["CHEBI:333"]


class TestParseEquationToPositionMap:
    """Tests for position-indexed equation parsing."""

    def test_simple_mapping(self):
        """Test simple position mapping."""
        left, right = parse_equation_to_position_map(
            "CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938"
        )
        assert left == {0: "CHEBI:16459", 1: "CHEBI:15377"}
        assert right == {0: "CHEBI:31011", 1: "CHEBI:28938"}

    def test_with_stoichiometry(self):
        """Test position mapping with stoichiometry expansion."""
        left, right = parse_equation_to_position_map(
            "2 CHEBI:15378 + CHEBI:30616 = CHEBI:456216"
        )
        # 2 CHEBI:15378 expands to positions 0, 1
        assert left[0] == "CHEBI:15378"
        assert left[1] == "CHEBI:15378"
        assert left[2] == "CHEBI:30616"
        assert right == {0: "CHEBI:456216"}

    def test_empty_equation(self):
        """Test empty equation returns empty maps."""
        left, right = parse_equation_to_position_map("")
        assert left == {}
        assert right == {}


class TestEquationChEBILookup:
    """Tests for EquationChEBILookup class."""

    def test_empty_lookup(self):
        """Test empty lookup initialization."""
        lookup = EquationChEBILookup()
        assert len(lookup) == 0
        assert lookup.get_equation("10000") is None
        assert lookup.get_chebi_ids("10000") == ([], [])

    def test_manual_equation_addition(self):
        """Test manually adding equations."""
        lookup = EquationChEBILookup()
        lookup._equations["10000"] = "CHEBI:16459 + CHEBI:15377 = CHEBI:31011"

        assert len(lookup) == 1
        assert lookup.get_equation("10000") == "CHEBI:16459 + CHEBI:15377 = CHEBI:31011"

        left, right = lookup.get_chebi_ids("10000")
        assert left == ["CHEBI:16459", "CHEBI:15377"]
        assert right == ["CHEBI:31011"]

    def test_caching(self):
        """Test that parsed results are cached."""
        lookup = EquationChEBILookup()
        lookup._equations["10000"] = "CHEBI:111 = CHEBI:222"

        # First call - should parse
        result1 = lookup.get_chebi_ids("10000")
        assert "10000" in lookup._parsed_cache

        # Second call - should use cache
        result2 = lookup.get_chebi_ids("10000")
        assert result1 == result2

    def test_load_from_file(self):
        """Test loading from gzipped file."""
        lookup = EquationChEBILookup()

        # Create a temporary gzipped file
        content = b"""ENTRY       RHEA:10000
DEFINITION  test reaction 1
EQUATION    CHEBI:16459 + CHEBI:15377 = CHEBI:31011
///
ENTRY       RHEA:10004
DEFINITION  test reaction 2
EQUATION    CHEBI:17484 = CHEBI:16017
///
ENTRY       RHEA:10008
DEFINITION  test reaction 3
EQUATION    CHEBI:35924 + CHEBI:29950 = CHEBI:50058
///
"""
        with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as f:
            with gzip.open(f.name, "wb") as gz:
                gz.write(content)
            temp_path = f.name

        try:
            count = lookup.load_from_file(temp_path)
            assert count == 3
            assert len(lookup) == 3

            # Check equations were loaded correctly
            assert lookup.get_equation("10000") == "CHEBI:16459 + CHEBI:15377 = CHEBI:31011"
            assert lookup.get_equation("10004") == "CHEBI:17484 = CHEBI:16017"
            assert lookup.get_equation("10008") == "CHEBI:35924 + CHEBI:29950 = CHEBI:50058"

            # Check parsing works
            left, right = lookup.get_chebi_ids("10000")
            assert left == ["CHEBI:16459", "CHEBI:15377"]
            assert right == ["CHEBI:31011"]

        finally:
            Path(temp_path).unlink()

    def test_nonexistent_rhea_id(self):
        """Test looking up non-existent RHEA ID."""
        lookup = EquationChEBILookup()
        lookup._equations["10000"] = "CHEBI:111 = CHEBI:222"

        assert lookup.get_equation("99999") is None
        assert lookup.get_chebi_ids("99999") == ([], [])


class TestCommonMoleculeChebiMapping:
    """Tests for the COMMON_MOLECULE_CHEBI dictionary completeness."""

    def test_essential_cofactors_present(self):
        """Ensure all essential cofactors are mapped."""
        essential = [
            "atp", "adp", "amp",
            "gtp", "gdp", "gmp",
            "nad(+)", "nadh", "nadp(+)", "nadph",
            "fad", "fadh2",
            "coa",
        ]
        for name in essential:
            assert name in COMMON_MOLECULE_CHEBI, f"Missing essential cofactor: {name}"

    def test_common_metabolites_present(self):
        """Ensure common metabolites are mapped."""
        metabolites = [
            "water", "phosphate", "h(+)",
            "o2", "co2", "h2",
            "pyruvate", "acetate",
        ]
        for name in metabolites:
            assert name in COMMON_MOLECULE_CHEBI, f"Missing metabolite: {name}"

    def test_all_chebi_ids_properly_formatted(self):
        """Ensure all ChEBI IDs have proper format."""
        for name, chebi_id in COMMON_MOLECULE_CHEBI.items():
            assert chebi_id.startswith("CHEBI:"), f"Invalid ChEBI format for {name}: {chebi_id}"
            # Check the number part is actually a number
            number_part = chebi_id[6:]
            assert number_part.isdigit(), f"Invalid ChEBI number for {name}: {chebi_id}"

    def test_no_duplicate_chebi_ids_for_unrelated_molecules(self):
        """Check that different molecule classes don't share ChEBI IDs (except aliases)."""
        # Group by ChEBI ID
        by_chebi: dict = {}
        for name, chebi_id in COMMON_MOLECULE_CHEBI.items():
            if chebi_id not in by_chebi:
                by_chebi[chebi_id] = []
            by_chebi[chebi_id].append(name)

        # Each ChEBI ID should only map to related names (aliases of same molecule)
        # This is informational - just check that we're aware of duplicates
        for chebi_id, names in by_chebi.items():
            if len(names) > 1:
                # All names should be related (this is a soft check)
                # Just make sure we don't have obviously unrelated molecules
                pass  # Could add specific checks if needed
