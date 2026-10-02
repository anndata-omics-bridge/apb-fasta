"""Peptide-to-protein matching and feature-aligned summaries."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import polars as pl
from prozor.matching.annotation import annotate_peptides_streaming

from apb_fasta.calculation.results import PeptideCoverage, PeptideLevelMatch


@dataclass(frozen=True, slots=True)
class PeptideLevelInput:
    """The columns needed to validate one peptide-derived feature level."""

    frame: pl.DataFrame
    sequence_column: str
    accession_column: str | None


@dataclass(frozen=True, slots=True)
class _ProteinRecord:
    id: str
    sequence: str


def match_peptide_levels(
    levels: Mapping[str, PeptideLevelInput],
    proteins: pl.DataFrame,
    /,
    *,
    backend: str,
    il_equivalent: bool,
    protein_group_separator: str,
) -> dict[str, PeptideLevelMatch]:
    """Match distinct normalized peptides once and summarize every source feature."""
    if not levels:
        return {}
    prepared = {name: _prepare_features(level, il_equivalent) for name, level in levels.items()}
    peptides = (
        pl.concat([frame.select("peptide") for frame in prepared.values()])
        .drop_nulls()
        .unique(maintain_order=True)
    )
    records = tuple(
        _ProteinRecord(
            id=str(index),
            sequence=_protein_sequence(value, il_equivalent),
        )
        for index, value in enumerate(proteins.get_column("sequence"))
    )
    matched = annotate_peptides_streaming(
        peptides.get_column("peptide").to_list(), records, backend=backend
    )
    occurrences = pl.DataFrame(
        [(match.peptide, int(match.protein_id)) for match in matched],
        schema={"peptide": pl.String, "record": pl.UInt32},
        orient="row",
    )
    site_counts = occurrences.group_by("peptide").agg(
        pl.len().cast(pl.UInt64).alias("fasta_match_site_count")
    )
    columns = set(proteins.columns)
    indexed = proteins.select(
        id=_text_column(proteins, "id"),
        organism=_text_column(
            proteins, "organism_mnemonic" if "organism_mnemonic" in columns else None
        ),
        is_contaminant=pl.col("is_contaminant") if "is_contaminant" in columns else pl.lit(False),
    ).with_row_index("record")
    distinct = occurrences.unique(maintain_order=True).join(
        indexed, on="record", how="left", maintain_order="left"
    )
    if distinct.get_column("id").null_count():
        raise ValueError("protein frame column 'id' contains a non-text value")
    matches = distinct.group_by("peptide", maintain_order=True).agg(
        pl.col("record").alias("matched_records"),
        pl.len().cast(pl.UInt64).alias("fasta_matching_protein_count"),
        pl.col("id").str.join(";").alias("fasta_matching_protein_ids"),
        fasta_matching_organisms=pl.col("organism").drop_nulls().unique().sort().str.join(";"),
        fasta_matches_contaminant=pl.col("is_contaminant").any(),
    )
    matches = matches.join(site_counts, on="peptide", how="left")
    assignments = _reported_assignments(prepared, proteins, protein_group_separator)
    return {
        name: _level_match(frame, matches, assignments, levels[name].accession_column is not None)
        for name, frame in prepared.items()
    }


def _text_column(frame: pl.DataFrame, column: str | None) -> pl.Expr:
    if column is None:
        return pl.lit(None, dtype=pl.String)
    dtype = frame.schema[column]
    if isinstance(dtype, (pl.String, pl.Categorical, pl.Enum)):
        return pl.col(column).cast(pl.String)
    if dtype == pl.Object:
        # Only foreign mixed-object columns need scalar type checking.
        return pl.col(column).map_elements(_optional_text, return_dtype=pl.String)
    return pl.lit(None, dtype=pl.String)


def _optional_text(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _prepare_features(level: PeptideLevelInput, il_equivalent: bool) -> pl.DataFrame:
    sequence = _text_column(level.frame, level.sequence_column).str.strip_chars().str.to_uppercase()
    if il_equivalent:
        sequence = sequence.str.replace_all("I", "L", literal=True)
    # with_columns keeps the source height even when both expressions are scalar nulls.
    return level.frame.with_columns(
        sequence.alias("peptide"),
        _text_column(level.frame, level.accession_column).alias("assignment"),
    ).select(
        pl.when(pl.col("peptide") != "").then(pl.col("peptide")).alias("peptide"),
        "assignment",
    )


def _reported_assignments(
    prepared: Mapping[str, pl.DataFrame], proteins: pl.DataFrame, separator: str
) -> pl.DataFrame:
    aliases = pl.concat(
        [
            proteins.with_columns(_text_column(proteins, column).alias("member"))
            .select("member")
            .with_row_index("record")
            for column in ("id", "accession")
            if column in proteins.columns
        ]
    )
    aliases = (
        aliases.filter(pl.col("member").is_not_null() & (pl.col("member") != ""))
        .unique()
        .group_by("member")
        .agg("record")
    )
    assignments = (
        pl.concat([frame.select("assignment") for frame in prepared.values()]).drop_nulls().unique()
    )
    return (
        assignments.with_columns(
            pl.col("assignment")
            .str.split(separator)
            .list.eval(pl.element().str.strip_chars())
            .alias("member")
        )
        .explode("member", empty_as_null=True)
        .filter(pl.col("member") != "")
        .join(aliases, on="member", how="left")
        .group_by("assignment")
        .agg(
            pl.len().cast(pl.UInt64).alias("reported_member_count"),
            pl.col("record")
            .is_not_null()
            .sum()
            .cast(pl.UInt64)
            .alias("reported_members_in_fasta_count"),
            pl.col("record")
            .explode(empty_as_null=True)
            .drop_nulls()
            .unique()
            .alias("reported_records"),
        )
    )


def _level_match(
    features: pl.DataFrame,
    matches: pl.DataFrame,
    assignments: pl.DataFrame,
    has_assignment: bool,
) -> PeptideLevelMatch:
    joined = features.join(matches, on="peptide", how="left", maintain_order="left").join(
        assignments, on="assignment", how="left", maintain_order="left"
    )
    summary = joined.select(
        (pl.col("fasta_matching_protein_count").fill_null(0) > 0).alias("peptide_in_fasta"),
        pl.col("fasta_match_site_count").fill_null(0),
        pl.col("fasta_matching_protein_count").fill_null(0),
        pl.col("fasta_matching_protein_ids").fill_null(""),
        pl.col("fasta_matching_organisms").fill_null(""),
        pl.col("fasta_matches_contaminant").fill_null(False),
        pl.col("reported_member_count").fill_null(0),
        pl.col("reported_members_in_fasta_count").fill_null(0),
        (
            pl.col("matched_records")
            .fill_null([])
            .list.set_intersection(pl.col("reported_records").fill_null([]))
            .list.len()
            > 0
        ).alias("peptide_in_reported_protein"),
    )
    if not has_assignment:
        summary = summary.with_columns(
            pl.lit(None, dtype=pl.UInt64).alias("reported_member_count"),
            pl.lit(None, dtype=pl.UInt64).alias("reported_members_in_fasta_count"),
            pl.lit(None, dtype=pl.Boolean).alias("peptide_in_reported_protein"),
        )
    unique = joined.select("peptide", "fasta_match_site_count").drop_nulls("peptide").unique()
    matched_count = int(summary.get_column("peptide_in_fasta").sum() or 0)
    return PeptideLevelMatch(
        summary=summary,
        coverage=PeptideCoverage(
            feature_count=features.height,
            unique_sequence_count=unique.height,
            matched_feature_count=matched_count,
            unmatched_feature_count=features.height - matched_count,
            match_site_count=int(unique.get_column("fasta_match_site_count").sum() or 0),
        ),
    )


def _protein_sequence(value: object, il_equivalent: bool) -> str:
    if not isinstance(value, str):
        raise ValueError("protein frame column 'sequence' contains a non-text value")
    return value.replace("I", "L") if il_equivalent else value
