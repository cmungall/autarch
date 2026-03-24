"""Unit tests for strict GO term validator alignment checks."""

from typing import Any

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.reaction import ReactionClass
from autarch.validation.go_term_validator import GoTermValidator


class FakeGoAdapter:
    """Simple fake adapter for GO label lookup in tests."""

    def __init__(self, label_map: dict[str, str]):
        self._label_map = label_map

    def label(self, go_id: str) -> str | None:
        return self._label_map.get(go_id)


class FakeEcAdapter:
    """Simple fake adapter for EC label lookup in tests."""

    def __init__(self, label_map: dict[str, str]):
        self._label_map = label_map

    def label(self, ec_id: str) -> str | None:
        return self._label_map.get(ec_id)


def make_reaction_class(
    class_name: str,
    module_name: str,
    go_id: str,
    docstring: str,
) -> type[ReactionClass]:
    """Create a minimal concrete ReactionClass for validation tests."""

    def check_membership_impl(self: ReactionClass, reaction: Reaction) -> ClassificationResult:
        return ClassificationResult(is_member=False, explanation="test")

    attrs: dict[str, Any] = {
        "__module__": module_name,
        "__doc__": docstring,
        "GO_ID": go_id,
        "check_membership_impl": check_membership_impl,
    }
    return type(class_name, (ReactionClass,), attrs)


def test_strict_alignment_passes_when_all_fields_match() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.rna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for RNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "RNA polymerase activity"}))

    result = validator.validate_class(cls)

    assert result.matches
    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_detects_filename_class_mismatch() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.dna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for RNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "RNA polymerase activity"}))

    result = validator.validate_class(cls)

    assert not result.matches
    assert not result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_detects_go_label_class_mismatch() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.rna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for RNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "DNA polymerase activity"}))

    result = validator.validate_class(cls)

    assert not result.matches
    assert result.file_matches_class
    assert not result.go_matches_class
    assert result.description_matches_class


def test_detects_description_class_mismatch() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.rna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for DNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "RNA polymerase activity"}))

    result = validator.validate_class(cls)

    assert not result.matches
    assert result.file_matches_class
    assert result.go_matches_class
    assert not result.description_matches_class


def test_wrapped_first_docstring_paragraph_still_matches_class() -> None:
    cls = make_reaction_class(
        class_name="OxidoreductaseActingOnTheCHOHGroupOfDonorsOxygenAsAcceptor",
        module_name="synthetic.oxidoreductase_acting_on_the_ch_oh_group_of_donors_oxygen_as_acceptor",
        go_id="GO:TEST",
        docstring=(
            "oxidoreductase acting on the CH-OH group of donors oxygen as\n"
            "acceptor.\n\n"
            "Catalysis of an oxidation-reduction reaction in which the donor is a\n"
            "CH-OH group and molecular oxygen is the acceptor."
        ),
    )
    validator = GoTermValidator(
        adapter=FakeGoAdapter(
            {"GO:TEST": "oxidoreductase activity, acting on the CH-OH group of donors, oxygen as acceptor"}
        )
    )

    result = validator.validate_class(cls)

    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_multiline_title_stops_before_definition_sentence() -> None:
    cls = make_reaction_class(
        class_name="ATPHydrolysis",
        module_name="synthetic.atp_hydrolysis",
        go_id="GO:TEST",
        docstring=(
            "ATP hydrolysis\n"
            "Catalysis of the reaction: ATP + H2O = ADP + H+ + phosphate.\n"
            "Examples:\n"
            "- ATPase"
        ),
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "ATP hydrolysis activity"}))

    result = validator.validate_class(cls)

    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_numeric_word_and_digit_forms_are_treated_equivalently() -> None:
    cls = make_reaction_class(
        class_name="TwoOxoglutarateDependentDioxygenase",
        module_name="synthetic.two_oxoglutarate_dependent_dioxygenase",
        go_id="GO:TEST",
        docstring="2-oxoglutarate-dependent dioxygenase.",
    )
    validator = GoTermValidator(
        adapter=FakeGoAdapter({"GO:TEST": "2-oxoglutarate-dependent dioxygenase activity"})
    )

    result = validator.validate_class(cls)

    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_acronym_suffixes_are_split_before_trailing_words() -> None:
    cls = make_reaction_class(
        class_name="ThreeBetaHydroxysteroid3DehydrogenaseNADPActivity",
        module_name="synthetic.three_beta_hydroxysteroid_3_dehydrogenase_nad_p_activity",
        go_id="GO:TEST",
        docstring="3-beta-hydroxysteroid 3-dehydrogenase (NADP+) activity.",
    )
    validator = GoTermValidator(
        adapter=FakeGoAdapter(
            {"GO:TEST": "3-beta-hydroxysteroid 3-dehydrogenase (NADP+) activity"}
        )
    )

    result = validator.validate_class(cls)

    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_teen_numbers_normalize_against_numeric_tokens() -> None:
    cls = make_reaction_class(
        class_name="FifteenOxoprostaglandin13ReductaseNADOrNADPActivity",
        module_name="synthetic.fifteen_oxoprostaglandin_13_reductase_nad_or_nadp_activity",
        go_id="GO:TEST",
        docstring="15-oxoprostaglandin 13-reductase [NAD(P)+] activity.",
    )
    validator = GoTermValidator(
        adapter=FakeGoAdapter(
            {"GO:TEST": "15-oxoprostaglandin 13-reductase [NAD(P)+] activity"}
        )
    )

    result = validator.validate_class(cls)

    assert result.file_matches_class
    assert result.go_matches_class
    assert result.description_matches_class


def test_ec_internal_labels_are_composed_with_parent_context() -> None:
    class TranslocationOfInorganicAnionsAndTheirChelatesLinkedToTheHydrolysisOfANucleosideTriphosphate(
        ReactionClass
    ):
        """translocation of inorganic anions and their chelates linked to the hydrolysis of a nucleoside triphosphate."""

        __module__ = (
            "synthetic.translocation_of_inorganic_anions_and_their_chelates_"
            "linked_to_the_hydrolysis_of_a_nucleoside_triphosphate"
        )
        GO_ID = None
        EC_NUMBER_PREFIX = "7.3.2.-"

        def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
            return ClassificationResult(is_member=False, explanation="test")

    validator = GoTermValidator(
        adapter=FakeGoAdapter({}),
        ec_adapter=FakeEcAdapter(
            {
                "EC:7.3": "Catalysing the translocation of inorganic anions and their chelates",
                "EC:7.3.2": "Linked to the hydrolysis of a nucleoside triphosphate",
            }
        ),
    )

    result = validator.validate_class(
        TranslocationOfInorganicAnionsAndTheirChelatesLinkedToTheHydrolysisOfANucleosideTriphosphate
    )

    assert result.go_matches_class
    assert result.matches


def test_ec_with_labels_are_composed_with_parent_context() -> None:
    class ActingOnOtherNitrogenousCompoundsAsDonorsWithACytochromeAsAcceptor(
        ReactionClass
    ):
        """acting on other nitrogenous compounds as donors with a cytochrome as acceptor."""

        __module__ = (
            "synthetic.acting_on_other_nitrogenous_compounds_as_donors_"
            "with_a_cytochrome_as_acceptor"
        )
        GO_ID = None
        EC_NUMBER_PREFIX = "1.7.2.-"

        def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
            return ClassificationResult(is_member=False, explanation="test")

    validator = GoTermValidator(
        adapter=FakeGoAdapter({}),
        ec_adapter=FakeEcAdapter(
            {
                "EC:1.7": "Acting on other nitrogenous compounds as donors",
                "EC:1.7.2": "With a cytochrome as acceptor",
            }
        ),
    )

    result = validator.validate_class(
        ActingOnOtherNitrogenousCompoundsAsDonorsWithACytochromeAsAcceptor
    )

    assert result.go_matches_class
    assert result.matches


def test_detects_broad_ec_prefix_against_exact_go_ec() -> None:
    cls = make_reaction_class(
        class_name="AlcoholDehydrogenaseNAD",
        module_name="synthetic.alcohol_dehydrogenase_nad",
        go_id="GO:TEST",
        docstring="Classifier for alcohol dehydrogenase NAD reactions.",
    )
    cls.EC_NUMBER_PREFIX = "1.1.1.-"
    validator = GoTermValidator(
        adapter=FakeGoAdapter({"GO:TEST": "alcohol dehydrogenase NAD activity"}),
        go_term_map={
            "GO:TEST": {
                "label": "alcohol dehydrogenase NAD activity",
                "definition": "Catalysis of the reaction: an alcohol + NAD+ = an aldehyde or ketone + NADH + H+.",
                "namespace": "molecular_function",
                "ec_numbers": ["1.1.1.1"],
                "ancestors": ["GO:TEST"],
            }
        },
    )

    result = validator.validate_class(cls)

    assert not result.matches
    assert not result.ec_matches_go
    assert "broader than exact GO EC" in result.explanation


def test_detects_mapping_driven_source_exclusion_anti_pattern() -> None:
    class AlcoholDehydrogenaseNAD(ReactionClass):
        """Classifier for alcohol dehydrogenase NAD reactions."""

        GO_ID = "GO:TEST"
        EC_NUMBER_PREFIX = "1.1.1.1"
        NON_4022_LABEL_EXCLUSIONS = ("xylitol",)

        def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
            return ClassificationResult(is_member=False, explanation="test")

    validator = GoTermValidator(
        adapter=FakeGoAdapter({"GO:TEST": "alcohol dehydrogenase NAD activity"}),
        go_term_map={
            "GO:TEST": {
                "label": "alcohol dehydrogenase NAD activity",
                "definition": "Catalysis of the reaction: an alcohol + NAD+ = an aldehyde or ketone + NADH + H+.",
                "namespace": "molecular_function",
                "ec_numbers": ["1.1.1.1"],
                "ancestors": ["GO:TEST"],
            }
        },
    )

    result = validator.validate_class(AlcoholDehydrogenaseNAD)

    assert not result.matches
    assert not result.source_matches_go
    assert "mapping-driven logic" in result.explanation


def test_detects_definition_conflict_with_ketone_context_gate() -> None:
    class AlcoholDehydrogenaseNAD(ReactionClass):
        """Classifier for alcohol dehydrogenase NAD reactions."""

        GO_ID = "GO:TEST"
        EC_NUMBER_PREFIX = "1.1.1.1"

        def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
            sugar_polyol_context = True
            if sugar_polyol_context:
                return ClassificationResult(is_member=False, explanation="test")
            return ClassificationResult(is_member=False, explanation="test")

    validator = GoTermValidator(
        adapter=FakeGoAdapter({"GO:TEST": "alcohol dehydrogenase NAD activity"}),
        go_term_map={
            "GO:TEST": {
                "label": "alcohol dehydrogenase NAD activity",
                "definition": "Catalysis of the reaction: an alcohol + NAD+ = an aldehyde or ketone + NADH + H+.",
                "namespace": "molecular_function",
                "ec_numbers": ["1.1.1.1"],
                "ancestors": ["GO:TEST"],
            }
        },
    )

    result = validator.validate_class(AlcoholDehydrogenaseNAD)

    assert not result.matches
    assert not result.definition_matches_go
    assert "aldehyde/ketone" in result.explanation


def test_detects_text_mining_source_anti_pattern() -> None:
    class AlcoholDehydrogenaseNAD(ReactionClass):
        """Classifier for alcohol dehydrogenase NAD reactions."""

        GO_ID = "GO:TEST"
        EC_NUMBER_PREFIX = "1.1.1.1"

        def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
            label_lower = reaction.label.lower() if reaction.label else ""
            if "alcohol" in label_lower:
                return ClassificationResult(is_member=True, explanation="text heuristic")
            return ClassificationResult(is_member=False, explanation="test")

    validator = GoTermValidator(
        adapter=FakeGoAdapter({"GO:TEST": "alcohol dehydrogenase NAD activity"}),
        go_term_map={
            "GO:TEST": {
                "label": "alcohol dehydrogenase NAD activity",
                "definition": "Catalysis of the reaction: an alcohol + NAD+ = an aldehyde or ketone + NADH + H+.",
                "namespace": "molecular_function",
                "ec_numbers": ["1.1.1.1"],
                "ancestors": ["GO:TEST"],
            }
        },
    )

    result = validator.validate_class(AlcoholDehydrogenaseNAD)

    assert not result.matches
    assert not result.source_matches_go
    assert "reaction.label text parsing" in result.explanation


def test_compare_against_baseline_ignores_known_failures() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.dna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for DNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "RNA polymerase activity"}))

    result = validator.validate_class(cls)
    baseline = {result.class_name: result.failure_flags}

    unexpected, resolved = validator.compare_against_baseline([result], baseline)

    assert unexpected == []
    assert resolved == []


def test_compare_against_baseline_detects_new_failures() -> None:
    cls = make_reaction_class(
        class_name="RNAPolymerase",
        module_name="synthetic.dna_polymerase",
        go_id="GO:TEST",
        docstring="Classifier for DNA polymerase reactions.",
    )
    validator = GoTermValidator(adapter=FakeGoAdapter({"GO:TEST": "RNA polymerase activity"}))

    result = validator.validate_class(cls)
    baseline = {result.class_name: ["file_matches_class"]}

    unexpected, resolved = validator.compare_against_baseline([result], baseline)

    assert len(unexpected) == 1
    assert unexpected[0].class_name == "RNAPolymerase"
    assert resolved == []
