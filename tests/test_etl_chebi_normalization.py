"""Tests for ChEBI normalization integration in RHEA ETL."""

import gzip

import pytest

from autarch.datamodel import Participant, Reaction
from autarch.etl.chebi_normalization import EquationChEBILookup
from autarch.etl.rhea_tsv_etl import RheaTSVETL


class TestNormalizeParticipantChebi:
    """Tests for normalize_participant_chebi method."""

    @pytest.fixture
    def etl(self):
        """Create ETL instance."""
        return RheaTSVETL()

    def test_keeps_existing_chebi_id(self, etl):
        """Existing ChEBI ID should not be overwritten."""
        p = Participant(name="test", chebi_id="CHEBI:12345")
        result = etl.normalize_participant_chebi(p)
        assert result.chebi_id == "CHEBI:12345"

    def test_uses_equation_chebi_id_by_position(self, etl):
        """ChEBI ID should come from EQUATION by position."""
        p = Participant(name="unknown")
        equation_ids = ["CHEBI:111", "CHEBI:222", "CHEBI:333"]

        result = etl.normalize_participant_chebi(p, equation_ids, position=0)
        assert result.chebi_id == "CHEBI:111"

        p2 = Participant(name="unknown2")
        result2 = etl.normalize_participant_chebi(p2, equation_ids, position=2)
        assert result2.chebi_id == "CHEBI:333"

    def test_uses_name_normalization_as_fallback(self, etl):
        """Name-based normalization should work when other methods fail."""
        p = Participant(name="H(+)")
        result = etl.normalize_participant_chebi(p)
        assert result.chebi_id == "CHEBI:15378"

        p2 = Participant(name="water")
        result2 = etl.normalize_participant_chebi(p2)
        assert result2.chebi_id == "CHEBI:15377"

        p3 = Participant(name="ATP")
        result3 = etl.normalize_participant_chebi(p3)
        assert result3.chebi_id == "CHEBI:30616"

    def test_equation_takes_priority_over_name(self, etl):
        """EQUATION ChEBI should be used even if name matches."""
        # Participant named "water" but EQUATION says different ChEBI
        p = Participant(name="water")
        equation_ids = ["CHEBI:99999"]

        result = etl.normalize_participant_chebi(p, equation_ids, position=0)
        assert result.chebi_id == "CHEBI:99999"  # EQUATION wins

    def test_position_out_of_range_falls_back_to_name(self, etl):
        """If position exceeds EQUATION list, fall back to name."""
        p = Participant(name="H(+)")
        equation_ids = ["CHEBI:111"]  # Only one ID

        result = etl.normalize_participant_chebi(p, equation_ids, position=5)
        # Position 5 is out of range, so should use name-based
        assert result.chebi_id == "CHEBI:15378"

    def test_returns_unchanged_if_no_match(self, etl):
        """Participant without matches should be unchanged."""
        p = Participant(name="some_unknown_molecule_xyz")
        result = etl.normalize_participant_chebi(p)
        assert result.chebi_id is None

    def test_empty_name_no_crash(self, etl):
        """Empty name should not crash."""
        p = Participant(name="")
        result = etl.normalize_participant_chebi(p)
        assert result.chebi_id is None

        p2 = Participant(name=None)
        result2 = etl.normalize_participant_chebi(p2)
        assert result2.chebi_id is None


class TestNormalizeReactionChebi:
    """Tests for normalize_reaction_chebi method."""

    @pytest.fixture
    def etl(self):
        """Create ETL instance."""
        return RheaTSVETL()

    @pytest.fixture
    def sample_reaction(self):
        """Create a sample reaction with missing ChEBI IDs."""
        return Reaction(
            left_participants=[
                Participant(name="ATP", chebi_id=None),
                Participant(name="H2O", chebi_id=None),
            ],
            right_participants=[
                Participant(name="ADP", chebi_id=None),
                Participant(name="phosphate", chebi_id=None),
                Participant(name="H(+)", chebi_id=None),
            ],
        )

    @pytest.fixture
    def equation_lookup(self):
        """Create a mock EQUATION lookup."""
        lookup = EquationChEBILookup()
        lookup._equations["10000"] = (
            "CHEBI:30616 + CHEBI:15377 = CHEBI:456216 + CHEBI:43474 + CHEBI:15378"
        )
        return lookup

    def test_normalizes_all_participants(self, etl, sample_reaction):
        """All participants with recognizable names should get ChEBI IDs."""
        etl.normalize_reaction_chebi(sample_reaction, "10000")

        # Check all participants got normalized via name lookup
        assert sample_reaction.left_participants[0].chebi_id == "CHEBI:30616"  # ATP
        assert sample_reaction.left_participants[1].chebi_id == "CHEBI:15377"  # H2O
        assert sample_reaction.right_participants[0].chebi_id == "CHEBI:456216"  # ADP
        assert sample_reaction.right_participants[1].chebi_id == "CHEBI:43474"  # phosphate
        assert sample_reaction.right_participants[2].chebi_id == "CHEBI:15378"  # H(+)

    def test_uses_equation_lookup(self, etl, equation_lookup):
        """EQUATION lookup should be used when available."""
        reaction = Reaction(
            left_participants=[
                Participant(name="unknown1", chebi_id=None),
                Participant(name="unknown2", chebi_id=None),
            ],
            right_participants=[
                Participant(name="unknown3", chebi_id=None),
            ],
        )

        etl.normalize_reaction_chebi(reaction, "10000", equation_lookup)

        # ChEBI IDs should come from EQUATION by position
        assert reaction.left_participants[0].chebi_id == "CHEBI:30616"
        assert reaction.left_participants[1].chebi_id == "CHEBI:15377"
        assert reaction.right_participants[0].chebi_id == "CHEBI:456216"

    def test_preserves_existing_chebi_ids(self, etl, equation_lookup):
        """Existing ChEBI IDs should not be overwritten."""
        reaction = Reaction(
            left_participants=[
                Participant(name="test", chebi_id="CHEBI:99999"),  # Already has ID
            ],
            right_participants=[
                Participant(name="unknown", chebi_id=None),
            ],
        )

        etl.normalize_reaction_chebi(reaction, "10000", equation_lookup)

        # Existing ID preserved
        assert reaction.left_participants[0].chebi_id == "CHEBI:99999"
        # Unknown gets from EQUATION
        assert reaction.right_participants[0].chebi_id == "CHEBI:456216"


class TestLoadEquationChEBILookup:
    """Tests for load_equation_chebi_lookup method."""

    @pytest.fixture
    def etl_with_temp_cache(self, tmp_path):
        """Create ETL with temporary cache directory."""
        etl = RheaTSVETL(cache_dir=tmp_path / "rhea_tsv")
        etl.cache_dir.mkdir(parents=True, exist_ok=True)

        # Create a minimal rhea-reactions.txt.gz
        content = b"""ENTRY       RHEA:10000
DEFINITION  pentanamide + H2O = pentanoate + NH4(+)
EQUATION    CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938
///
ENTRY       RHEA:10004
DEFINITION  benzyl isothiocyanate = benzyl thiocyanate
EQUATION    CHEBI:17484 = CHEBI:16017
///
"""
        gz_path = etl.cache_dir / "rhea-reactions.txt.gz"
        with gzip.open(gz_path, "wb") as f:
            f.write(content)

        return etl

    def test_loads_equations(self, etl_with_temp_cache):
        """Should load EQUATION data from file."""
        lookup = etl_with_temp_cache.load_equation_chebi_lookup()

        assert len(lookup) == 2
        assert lookup.get_equation("10000") == "CHEBI:16459 + CHEBI:15377 = CHEBI:31011 + CHEBI:28938"
        assert lookup.get_equation("10004") == "CHEBI:17484 = CHEBI:16017"

    def test_parsed_chebi_ids(self, etl_with_temp_cache):
        """Should correctly parse ChEBI IDs from equations."""
        lookup = etl_with_temp_cache.load_equation_chebi_lookup()

        left, right = lookup.get_chebi_ids("10000")
        assert left == ["CHEBI:16459", "CHEBI:15377"]
        assert right == ["CHEBI:31011", "CHEBI:28938"]


class TestCommonMoleculeNormalization:
    """Integration tests for common molecule normalization."""

    @pytest.fixture
    def etl(self):
        """Create ETL instance."""
        return RheaTSVETL()

    @pytest.mark.parametrize(
        "name,expected_chebi",
        [
            # Protons
            ("H(+)", "CHEBI:15378"),
            ("H+", "CHEBI:15378"),
            ("hydron", "CHEBI:15378"),
            # Water
            ("H2O", "CHEBI:15377"),
            ("water", "CHEBI:15377"),
            # Phosphate
            ("phosphate", "CHEBI:43474"),
            ("Pi", "CHEBI:43474"),
            # ATP family
            ("ATP", "CHEBI:30616"),
            ("ADP", "CHEBI:456216"),
            ("AMP", "CHEBI:456215"),
            # NAD cofactors
            ("NAD(+)", "CHEBI:57540"),
            ("NADH", "CHEBI:57945"),
            ("NADP(+)", "CHEBI:58349"),
            ("NADPH", "CHEBI:57783"),
            # Gases
            ("O2", "CHEBI:15379"),
            ("CO2", "CHEBI:16526"),
            ("NH3", "CHEBI:16134"),
            # Common metabolites
            ("acetate", "CHEBI:30089"),
            ("pyruvate", "CHEBI:15361"),
        ],
    )
    def test_common_molecules_normalize(self, etl, name, expected_chebi):
        """Common molecules should normalize to correct ChEBI IDs."""
        p = Participant(name=name)
        result = etl.normalize_participant_chebi(p)
        assert result.chebi_id == expected_chebi, f"Failed for {name}"

    def test_case_insensitivity(self, etl):
        """Normalization should be case-insensitive."""
        for name in ["atp", "ATP", "Atp", "aTp"]:
            p = Participant(name=name)
            result = etl.normalize_participant_chebi(p)
            assert result.chebi_id == "CHEBI:30616", f"Failed for {name}"
