"""Peptide-to-protein matching and feature-aligned summaries."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise, product

import polars as pl
from prozor.api import annotate_peptides

from apb_fasta.calculation.results import PeptideCoverage, PeptideLevelMatch

# prolfquapp's tryptic peptides: cleaved after K or R unless P follows, 7 to 29 residues long.
_TRYPTIC_CLEAVAGE = re.compile(r"[KR](?!P|$)")
_TRYPTIC_LENGTHS = range(7, 30)
_LEADING_TEXT = {
    "reported_leading_id": "id",
    "reported_leading_accession": "accession",
    "reported_leading_description": "description",
    "reported_leading_gene_name": "gene_name",
}


@dataclass(frozen=True, slots=True)
class PeptideLevelInput:
    """The columns needed to validate one peptide-derived feature level."""

    frame: pl.DataFrame
    sequence_column: str
    accession_column: str | None
    decoy_column: str
    protein_assignment_column: str | None
    """The feature's protein group; its first member is the leading protein."""


def match_peptide_levels(
    levels: Mapping[str, PeptideLevelInput],
    proteins: pl.DataFrame,
    /,
    *,
    backend: str,
    protein_group_separator: str,
) -> dict[str, PeptideLevelMatch]:
    """Match distinct peptides exactly, then try other I/L spellings for the misses.

    Mass spectrometry cannot tell isoleucine from leucine, so a vendor may spell a peptide
    differently from its FASTA protein. Protein sequences are matched as read. A peptide
    with an exact match keeps exactly its proteins; its other I/L spellings are tried only
    when it has no exact match, or none in a reported protein the FASTA contains.
    """
    if not levels:
        return {}
    prepared = {
        name: _prepare_features(level, protein_group_separator) for name, level in levels.items()
    }
    peptides = (
        pl.concat([frame.select("peptide") for frame in prepared.values()])
        .drop_nulls()
        .unique(maintain_order=True)
        .get_column("peptide")
    )
    sequences = proteins.select(_text_column(proteins, "sequence")).to_series()
    if sequences.null_count():
        raise ValueError("protein frame column 'sequence' contains a non-text value")
    database = dict(zip(map(str, range(sequences.len())), sequences.to_list(), strict=True))
    exact = _occurrences(peptides.to_list(), database, backend)
    members = _member_records(proteins)
    assignments = _reported_assignments(prepared, members, protein_group_separator)
    leading = _leading_proteins(prepared, proteins, members)
    spellings = _il_spellings(_fallback_peptides(prepared, peptides, exact, assignments))
    il = (
        _occurrences(
            spellings.get_column("spelling").unique(maintain_order=True).to_list(),
            database,
            backend,
        )
        .rename({"peptide": "spelling"})
        .join(spellings, on="spelling")
        .select("peptide", "record")
    )
    occurrences = pl.concat(
        [
            exact.with_columns(il_only_match=pl.lit(value=False)),
            il.join(exact, on="peptide", how="anti").with_columns(il_only_match=pl.lit(value=True)),
        ]
    )
    site_counts = occurrences.group_by("peptide").agg(
        pl.len().cast(pl.UInt64).alias("fasta_match_site_count"),
        pl.col("il_only_match").first(),
    )
    il_records = il.unique().group_by("peptide").agg(pl.col("record").alias("il_records"))
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
    matches = matches.join(site_counts, on="peptide", how="left").join(
        il_records, on="peptide", how="left"
    )
    return {
        name: _level_match(
            frame, matches, assignments, leading, levels[name].accession_column is not None
        )
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


def _occurrences(peptides: list[str], database: dict[str, str], backend: str) -> pl.DataFrame:
    """One row per site where a peptide occurs in a protein sequence."""
    matched = annotate_peptides(peptides, database, backend=backend)
    return pl.DataFrame(
        [(match.peptide, int(match.protein_id)) for match in matched],
        schema={"peptide": pl.String, "record": pl.UInt32},
        orient="row",
    )


def _fallback_peptides(
    prepared: Mapping[str, pl.DataFrame],
    peptides: pl.Series,
    exact: pl.DataFrame,
    assignments: pl.DataFrame,
) -> pl.Series:
    """Peptides absent from the FASTA, or absent from every reported protein it contains."""
    exact_records = exact.unique().group_by("peptide").agg(pl.col("record").alias("exact_records"))
    outside_reported = (
        pl.concat([frame.select("peptide", "assignment") for frame in prepared.values()])
        .drop_nulls()
        .unique()
        .join(exact_records, on="peptide")
        .join(assignments, on="assignment")
        .filter(
            pl.col("reported_records").list.len() > 0,
            pl.col("exact_records").list.set_intersection("reported_records").list.len() == 0,
        )
        .get_column("peptide")
    )
    absent = (
        peptides.to_frame()
        .join(exact_records, on="peptide", how="anti", maintain_order="left")
        .get_column("peptide")
    )
    return pl.concat([absent, outside_reported]).unique(maintain_order=True)


def _il_spellings(peptides: pl.Series) -> pl.DataFrame:
    """Every other I/L spelling of each peptide, keyed by the peptide it stands for."""
    return pl.DataFrame(
        [
            (spelling, peptide)
            for peptide in peptides.to_list()
            for spelling in map(
                "".join,
                product(*(("I", "L") if residue in "IL" else residue for residue in peptide)),
            )
            if spelling != peptide
        ],
        schema={"spelling": pl.String, "peptide": pl.String},
        orient="row",
    )


def _prepare_features(level: PeptideLevelInput, separator: str) -> pl.DataFrame:
    sequence = _text_column(level.frame, level.sequence_column).str.strip_chars().str.to_uppercase()
    leading = (
        _text_column(level.frame, level.protein_assignment_column)
        .str.split(separator)
        .list.first()
        .str.strip_chars()
    )
    # with_columns keeps the source height even when the expressions are scalar nulls.
    return level.frame.with_columns(
        sequence.alias("peptide"),
        _text_column(level.frame, level.accession_column).alias("assignment"),
        leading.alias("leading_member"),
    ).select(
        pl.when(pl.col("peptide") != "").then(pl.col("peptide")).alias("peptide"),
        "assignment",
        pl.when(pl.col("leading_member") != "")
        .then(pl.col("leading_member"))
        .alias("leading_member"),
        pl.col(level.decoy_column).fill_null(value=False).alias("decoy"),
    )


def _member_records(proteins: pl.DataFrame) -> pl.DataFrame:
    """The protein records each reported member names, through their id or accession."""
    aliases = pl.concat(
        [
            proteins.with_columns(_text_column(proteins, column).alias("member"))
            .select("member")
            .with_row_index("record")
            for column in ("id", "accession")
            if column in proteins.columns
        ]
    )
    return (
        aliases.filter(pl.col("member").is_not_null() & (pl.col("member") != ""))
        .unique()
        .group_by("member")
        .agg("record")
    )


def _leading_proteins(
    prepared: Mapping[str, pl.DataFrame], proteins: pl.DataFrame, members: pl.DataFrame
) -> pl.DataFrame:
    """The first FASTA record each leading member names, as prolfquapp annotates a protein."""
    columns = set(proteins.columns)
    records = proteins.select(
        *(
            _text_column(proteins, column if column in columns else None).alias(name)
            for name, column in _LEADING_TEXT.items()
        ),
        _text_column(proteins, "sequence").alias("sequence"),
    ).with_row_index("record")
    return (
        pl.concat([frame.select("leading_member") for frame in prepared.values()])
        .drop_nulls()
        .unique()
        .join(members, left_on="leading_member", right_on="member")
        .select("leading_member", pl.col("record").list.min())
        .join(records, on="record")
        .select(
            "leading_member",
            *_LEADING_TEXT,
            pl.col("sequence")
            .str.len_chars()
            .cast(pl.UInt64)
            .alias("reported_leading_protein_length"),
            pl.col("sequence")
            .map_elements(_tryptic_peptides, return_dtype=pl.UInt64)
            .alias("reported_leading_tryptic_peptides"),
        )
    )


def _tryptic_peptides(sequence: str) -> int:
    ends = [match.end() for match in _TRYPTIC_CLEAVAGE.finditer(sequence.upper())]
    ends.append(len(sequence))
    return sum(end - start in _TRYPTIC_LENGTHS for start, end in pairwise([0, *ends]))


def _reported_assignments(
    prepared: Mapping[str, pl.DataFrame], aliases: pl.DataFrame, separator: str
) -> pl.DataFrame:
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
    leading: pl.DataFrame,
    has_assignment: bool,
) -> PeptideLevelMatch:
    joined = (
        features.join(matches, on="peptide", how="left", maintain_order="left")
        .join(assignments, on="assignment", how="left", maintain_order="left")
        .join(leading, on="leading_member", how="left", maintain_order="left")
    )
    reported = pl.col("reported_records").fill_null([])
    in_matched = pl.col("matched_records").fill_null([]).list.set_intersection(reported).list.len()
    in_il = pl.col("il_records").fill_null([]).list.set_intersection(reported).list.len()
    il_only_match = pl.col("il_only_match").fill_null(value=False)
    summary = joined.select(
        (pl.col("fasta_matching_protein_count").fill_null(0) > 0).alias("peptide_in_fasta"),
        pl.col("fasta_match_site_count").fill_null(0),
        pl.col("fasta_matching_protein_count").fill_null(0),
        pl.col("fasta_matching_protein_ids").fill_null(""),
        pl.col("fasta_matching_organisms").fill_null(""),
        pl.col("fasta_matches_contaminant").fill_null(False),
        pl.col("reported_member_count").fill_null(0),
        pl.col("reported_members_in_fasta_count").fill_null(0),
        ((in_matched > 0) | (in_il > 0)).alias("peptide_in_reported_protein"),
        (il_only_match | ((in_il > 0) & (in_matched == 0))).alias("fasta_il_only"),
        *leading.columns[1:],
    )
    if not has_assignment:
        summary = summary.with_columns(
            pl.lit(None, dtype=pl.UInt64).alias("reported_member_count"),
            pl.lit(None, dtype=pl.UInt64).alias("reported_members_in_fasta_count"),
            pl.lit(None, dtype=pl.Boolean).alias("peptide_in_reported_protein"),
        )
    unique = joined.select("peptide", "fasta_match_site_count").drop_nulls("peptide").unique()
    decoys = features.get_column("decoy")
    decoy_count = int(decoys.sum() or 0)
    matched_count = int(summary.get_column("peptide_in_fasta").filter(~decoys).sum() or 0)
    il_only_count = int(joined.select(il_only_match).to_series().filter(~decoys).sum() or 0)
    return PeptideLevelMatch(
        summary=summary,
        coverage=PeptideCoverage(
            feature_count=features.height,
            unique_sequence_count=unique.height,
            matched_feature_count=matched_count,
            unmatched_feature_count=features.height - decoy_count - matched_count,
            match_site_count=int(unique.get_column("fasta_match_site_count").sum() or 0),
            decoy_feature_count=decoy_count,
            il_only_matched_feature_count=il_only_count,
        ),
    )
