"""Native feature summaries preserve record identity, multiplicity and axis order."""

from __future__ import annotations

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from apb_fasta.calculation.matching import PeptideLevelInput, match_peptide_levels
from apb_fasta.calculation.results import PeptideCoverage, PeptideLevelMatch

_LEADING = [
    "reported_leading_id",
    "reported_leading_accession",
    "reported_leading_description",
    "reported_leading_gene_name",
    "reported_leading_protein_length",
    "reported_leading_tryptic_peptides",
]


def _proteins() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "id": ["p1", "p1", "p3"],
            "accession": ["a", "a", "p1"],
            "sequence": ["AAAIL", "AAA", "PILL"],
        }
    )


def _match(
    frame: pl.DataFrame,
    *,
    assignment: str | None = "assignment",
    leading: str | None = None,
    proteins: pl.DataFrame | None = None,
    backend: str = "auto",
    separator: str = ";",
) -> PeptideLevelMatch:
    if "decoy" not in frame.columns:
        frame = frame.with_columns(pl.lit(value=False).alias("decoy"))
    return match_peptide_levels(
        {"ion": PeptideLevelInput(frame, "peptide", assignment, "decoy", leading)},
        _proteins() if proteins is None else proteins,
        backend=backend,
        protein_group_separator=separator,
    )["ion"]


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_duplicate_records_sites_assignments_and_feature_order(backend: str) -> None:
    frame = pl.DataFrame(
        {
            "peptide": [" aa ", "AAA", None, "", "il", "pill", "ABSENT", "aa"],
            "assignment": ["a; a; unknown", "p1", None, "", "p1", "p3", "a", " ; "],
        }
    )
    result = _match(frame, backend=backend)
    expected = pl.DataFrame(
        {
            "peptide_in_fasta": [True, True, False, False, True, True, False, True],
            "fasta_match_site_count": [4, 2, 0, 0, 2, 1, 0, 4],
            "fasta_matching_protein_count": [2, 2, 0, 0, 2, 1, 0, 2],
            "fasta_matching_protein_ids": ["p1;p1", "p1;p1", "", "", "p1;p3", "p3", "", "p1;p1"],
            "fasta_matching_organisms": [""] * 8,
            "fasta_matches_contaminant": [False] * 8,
            "reported_member_count": [3, 1, 0, 0, 1, 1, 1, 0],
            "reported_members_in_fasta_count": [2, 1, 0, 0, 1, 1, 1, 0],
            "peptide_in_reported_protein": [True, True, False, False, True, True, False, False],
            "fasta_il_only": [False] * 8,
        },
        schema={
            "peptide_in_fasta": pl.Boolean,
            "fasta_match_site_count": pl.UInt64,
            "fasta_matching_protein_count": pl.UInt64,
            "fasta_matching_protein_ids": pl.String,
            "fasta_matching_organisms": pl.String,
            "fasta_matches_contaminant": pl.Boolean,
            "reported_member_count": pl.UInt64,
            "reported_members_in_fasta_count": pl.UInt64,
            "peptide_in_reported_protein": pl.Boolean,
            "fasta_il_only": pl.Boolean,
        },
    )
    assert result.summary.columns == [*expected.columns, *_LEADING]
    assert_frame_equal(result.summary.select(expected.columns), expected)
    assert result.coverage == PeptideCoverage(8, 5, 5, 3, 9, 0, 0)


def test_matched_organisms_and_contaminants_are_summarized_per_feature() -> None:
    proteins = pl.DataFrame(
        {
            "id": ["sp|Cont_P1|X_BOVIN", "sp|P1|X_HUMAN", "sp|P2|Y_YEAST"],
            "sequence": ["AAKK", "AAKR", "AAKR"],
            "organism_mnemonic": ["BOVIN", "HUMAN", "YEAST"],
            "is_contaminant": [True, False, False],
        }
    )
    result = _match(
        pl.DataFrame({"peptide": ["AAK", "AKR", "MM"]}), assignment=None, proteins=proteins
    )

    assert result.summary.get_column("fasta_matching_organisms").to_list() == [
        "BOVIN;HUMAN;YEAST",
        "HUMAN;YEAST",
        "",
    ]
    assert result.summary.get_column("fasta_matches_contaminant").to_list() == [True, False, False]


def test_absent_assignment_is_null_not_zero() -> None:
    result = _match(pl.DataFrame({"peptide": ["AA", None]}), assignment=None)
    for column in (
        "reported_member_count",
        "reported_members_in_fasta_count",
        "peptide_in_reported_protein",
    ):
        assert result.summary.get_column(column).null_count() == 2
    assert result.summary.get_column("peptide_in_fasta").to_list() == [True, False]


@pytest.mark.parametrize("size", [0, 2])
def test_empty_and_null_features_keep_height_and_schema(size: int) -> None:
    frame = pl.DataFrame({"peptide": [None] * size, "assignment": [None] * size})
    result = _match(frame)
    assert result.summary.height == size
    assert result.summary.schema["fasta_matching_protein_ids"] == pl.String
    assert result.summary.schema["fasta_match_site_count"] == pl.UInt64
    assert result.coverage == PeptideCoverage(size, 0, 0, size, 0, 0, 0)


def test_empty_database_and_missing_accession_column() -> None:
    frame = pl.DataFrame({"peptide": ["AA"], "assignment": ["p1"]})
    empty = _match(frame, proteins=_proteins().clear())
    assert empty.coverage == PeptideCoverage(1, 1, 0, 1, 0, 0, 0)
    assert empty.summary["reported_members_in_fasta_count"].to_list() == [0]
    no_accessions = _match(frame, proteins=_proteins().drop("accession"))
    assert no_accessions.summary["peptide_in_reported_protein"].to_list() == [True]


@pytest.mark.parametrize("dtype", [pl.String(), pl.Categorical(), pl.Enum(["aa", "a"])])
def test_string_like_inputs(dtype: pl.DataType) -> None:
    frame = pl.DataFrame(
        {"peptide": ["aa", None], "assignment": ["a", None]},
        schema_overrides={"peptide": dtype, "assignment": dtype},
    )
    assert _match(frame).summary["peptide_in_reported_protein"].to_list() == [True, False]


def test_mixed_object_and_nontext_inputs() -> None:
    frame = pl.DataFrame(
        {
            "peptide": pl.Series(["aa", 12, None], dtype=pl.Object),
            "assignment": pl.Series(["p1", 12, None], dtype=pl.Object),
        }
    )
    assert _match(frame).summary["peptide_in_reported_protein"].to_list() == [True, False, False]
    nontext = _match(pl.DataFrame({"peptide": [1, 2], "assignment": [1, 2]}))
    assert nontext.coverage == PeptideCoverage(2, 0, 0, 2, 0, 0, 0)


def test_decoys_count_apart_from_matched_and_unmatched_targets() -> None:
    frame = pl.DataFrame(
        {
            "peptide": ["AA", "GG", "GG", "AA"],
            "assignment": ["p1", None, None, None],
            "decoy": [False, False, True, True],
        }
    )

    coverage = _match(frame).coverage

    assert coverage.feature_count == 4
    assert (coverage.matched_feature_count, coverage.unmatched_feature_count) == (1, 1)
    assert coverage.decoy_feature_count == 2


def test_multilevel_coverage_counts_sequences_per_level() -> None:
    frame = pl.DataFrame(
        {"peptide": ["AA", "aa"], "assignment": ["a|a", "p1"], "decoy": [False, False]}
    )
    levels = {
        "ion": PeptideLevelInput(frame, "peptide", "assignment", "decoy", None),
        "peptide": PeptideLevelInput(frame.head(1), "peptide", "assignment", "decoy", None),
    }
    result = match_peptide_levels(levels, _proteins(), backend="auto", protein_group_separator="|")
    assert result["ion"].coverage == PeptideCoverage(2, 1, 2, 0, 4, 0, 0)
    assert result["peptide"].coverage == PeptideCoverage(1, 1, 1, 0, 4, 0, 0)
    assert result["ion"].summary["reported_member_count"].to_list() == [2, 1]


@pytest.mark.parametrize("column", ["id", "sequence"])
def test_nontext_protein_fields_are_rejected(column: str) -> None:
    proteins = _proteins().with_columns(pl.lit(12).alias(column))
    with pytest.raises(ValueError, match=f"column '{column}' contains a non-text value"):
        _match(pl.DataFrame({"peptide": ["AA"], "assignment": ["p1"]}), proteins=proteins)


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_other_il_spellings_fill_only_what_exact_matching_misses(backend: str) -> None:
    proteins = pl.DataFrame(
        {"id": ["P1", "P2", "P3"], "sequence": ["MKAILEK", "MKALLEK", "GGLIVR"]}
    )
    frame = pl.DataFrame(
        {
            "peptide": ["AILE", "GGILVR", "WWW", "ALLE", "ALLE", "GGILVR"],
            "assignment": ["P1", "P3", None, "P1", "P2", None],
            "decoy": [False, False, False, False, False, True],
        }
    )

    result = _match(frame, proteins=proteins, backend=backend)

    summary = result.summary
    # AILE occurs exactly in P1 and keeps it alone, although P2 has the spelling ALLE.
    assert summary["fasta_matching_protein_ids"].to_list() == ["P1", "P3", "", "P2", "P2", "P3"]
    assert summary["peptide_in_fasta"].to_list() == [True, True, False, True, True, True]
    # ALLE assigned to P1 occurs there only as AILE; its protein list stays exact.
    assert summary["peptide_in_reported_protein"].to_list() == [
        True,
        True,
        False,
        True,
        True,
        False,
    ]
    assert summary["fasta_il_only"].to_list() == [False, True, False, True, False, True]
    assert result.coverage == PeptideCoverage(6, 4, 4, 1, 3, 1, 1)


def test_leading_protein_is_the_first_reported_member_found_by_id_or_accession() -> None:
    proteins = pl.DataFrame(
        {
            "id": ["sp|P1|ONE_HUMAN", "sp|P2|TWO_HUMAN", "sp|P2|TWO_HUMAN"],
            "accession": ["P1", "P2", "P2"],
            "description": ["One OS=Homo sapiens GN=ONE", "Two", "Two again"],
            "gene_name": ["ONE", None, None],
            "sequence": ["MKGLPRAKSHGSTGWGKRKRNKPK", "MPEPTIDEKAACLLKR", "AAAA"],
        }
    )
    frame = pl.DataFrame(
        {
            "peptide": ["GLPR", "AACLLK", "AACLLK", "GLPR", "GLPR"],
            "group": [" P2 ; P1", "sp|P1|ONE_HUMAN", "unknown;P1", None, ""],
        }
    )

    summary = _match(frame, assignment=None, leading="group", proteins=proteins).summary

    expected = pl.DataFrame(
        {
            "reported_leading_id": ["sp|P2|TWO_HUMAN", "sp|P1|ONE_HUMAN", None, None, None],
            "reported_leading_accession": ["P2", "P1", None, None, None],
            "reported_leading_description": ["Two", "One OS=Homo sapiens GN=ONE", None, None, None],
            "reported_leading_gene_name": [None, "ONE", None, None, None],
            "reported_leading_protein_length": [16, 24, None, None, None],
            "reported_leading_tryptic_peptides": [1, 1, None, None, None],
        },
        schema_overrides={
            "reported_leading_gene_name": pl.String,
            "reported_leading_protein_length": pl.UInt64,
            "reported_leading_tryptic_peptides": pl.UInt64,
        },
    )
    assert_frame_equal(summary.select(_LEADING), expected)


@pytest.mark.parametrize(
    ("sequence", "count"),
    [
        # prolfquapp's nr_tryptic_peptides(min_length = 7, max_length = 30) on each sequence
        ("MKGLPRAKSHGSTGWGKRKRNKPK", 1),
        ("AAAAAAAKPAAAAAAARAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKAAAAAAAK", 2),
        ("kkkkaaaaaaar", 1),
        ("AAAAAA", 0),
    ],
)
def test_leading_protein_counts_tryptic_peptides_as_prolfquapp(sequence: str, count: int) -> None:
    proteins = pl.DataFrame({"id": ["P1"], "sequence": [sequence]})
    frame = pl.DataFrame({"peptide": ["AA"], "group": ["P1"]})

    summary = _match(frame, assignment=None, leading="group", proteins=proteins).summary

    assert summary["reported_leading_tryptic_peptides"].to_list() == [count]
    assert summary["reported_leading_accession"].to_list() == [None]
