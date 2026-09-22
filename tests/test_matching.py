"""Native feature summaries preserve record identity, multiplicity and axis order."""

from __future__ import annotations

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from apb_fasta.calculation.matching import PeptideLevelInput, match_peptide_levels
from apb_fasta.calculation.results import PeptideCoverage, PeptideLevelMatch


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
    proteins: pl.DataFrame | None = None,
    il_equivalent: bool = False,
    backend: str = "auto",
    separator: str = ";",
) -> PeptideLevelMatch:
    return match_peptide_levels(
        {"ion": PeptideLevelInput(frame, "peptide", assignment)},
        _proteins() if proteins is None else proteins,
        backend=backend,
        il_equivalent=il_equivalent,
        protein_group_separator=separator,
    )["ion"]


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
@pytest.mark.parametrize("il_equivalent", [False, True])
def test_duplicate_records_sites_assignments_and_feature_order(
    backend: str, il_equivalent: bool
) -> None:
    frame = pl.DataFrame(
        {
            "peptide": [" aa ", "AAA", None, "", "il", "pill", "ABSENT", "aa"],
            "assignment": ["a; a; unknown", "p1", None, "", "p1", "p3", "a", " ; "],
        }
    )
    result = _match(frame, backend=backend, il_equivalent=il_equivalent)
    expected = pl.DataFrame(
        {
            "peptide_in_fasta": [True, True, False, False, True, True, False, True],
            "fasta_match_site_count": [4, 2, 0, 0, 3 if il_equivalent else 2, 1, 0, 4],
            "fasta_matching_protein_count": [2, 2, 0, 0, 2, 1, 0, 2],
            "fasta_matching_protein_ids": ["p1;p1", "p1;p1", "", "", "p1;p3", "p3", "", "p1;p1"],
            "reported_member_count": [3, 1, 0, 0, 1, 1, 1, 0],
            "reported_members_in_fasta_count": [2, 1, 0, 0, 1, 1, 1, 0],
            "peptide_in_reported_protein": [True, True, False, False, True, True, False, False],
        },
        schema={
            "peptide_in_fasta": pl.Boolean,
            "fasta_match_site_count": pl.UInt64,
            "fasta_matching_protein_count": pl.UInt64,
            "fasta_matching_protein_ids": pl.String,
            "reported_member_count": pl.UInt64,
            "reported_members_in_fasta_count": pl.UInt64,
            "peptide_in_reported_protein": pl.Boolean,
        },
    )
    assert_frame_equal(result.summary, expected)
    assert result.coverage == PeptideCoverage(8, 5, 5, 3, 10 if il_equivalent else 9)


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
    assert result.coverage == PeptideCoverage(size, 0, 0, size, 0)


def test_empty_database_and_missing_accession_column() -> None:
    frame = pl.DataFrame({"peptide": ["AA"], "assignment": ["p1"]})
    empty = _match(frame, proteins=_proteins().clear())
    assert empty.coverage == PeptideCoverage(1, 1, 0, 1, 0)
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
    assert nontext.coverage == PeptideCoverage(2, 0, 0, 2, 0)


def test_multilevel_coverage_counts_sequences_per_level() -> None:
    frame = pl.DataFrame({"peptide": ["AA", "aa"], "assignment": ["a|a", "p1"]})
    levels = {
        "ion": PeptideLevelInput(frame, "peptide", "assignment"),
        "peptide": PeptideLevelInput(frame.head(1), "peptide", "assignment"),
    }
    result = match_peptide_levels(
        levels, _proteins(), backend="auto", il_equivalent=False, protein_group_separator="|"
    )
    assert result["ion"].coverage == PeptideCoverage(2, 1, 2, 0, 4)
    assert result["peptide"].coverage == PeptideCoverage(1, 1, 1, 0, 4)
    assert result["ion"].summary["reported_member_count"].to_list() == [2, 1]


@pytest.mark.parametrize("column", ["id", "sequence"])
def test_nontext_protein_fields_are_rejected(column: str) -> None:
    proteins = _proteins().with_columns(pl.lit(12).alias(column))
    with pytest.raises(ValueError, match=f"column '{column}' contains a non-text value"):
        _match(pl.DataFrame({"peptide": ["AA"], "assignment": ["p1"]}), proteins=proteins)
