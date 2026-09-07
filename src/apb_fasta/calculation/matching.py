"""Peptide-to-protein matching and feature-aligned summaries."""

from __future__ import annotations

from collections import defaultdict
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
    prepared = {
        name: tuple(
            _normalize_sequence(value, il_equivalent)
            for value in level.frame[level.sequence_column]
        )
        for name, level in levels.items()
    }
    peptides = tuple(
        dict.fromkeys(
            peptide for values in prepared.values() for peptide in values if peptide is not None
        )
    )
    protein_rows = proteins.to_dicts()
    records = tuple(
        _ProteinRecord(
            id=str(index),
            sequence=_equivalent_sequence(_required_text(row, "sequence"), il_equivalent),
        )
        for index, row in enumerate(protein_rows)
    )
    matched = annotate_peptides_streaming(peptides, records, backend=backend)
    occurrences: dict[str, list[int]] = defaultdict(list)
    site_counts: dict[str, int] = defaultdict(int)
    for match in matched:
        record_index = int(match.protein_id)
        occurrences[match.peptide].append(record_index)
        site_counts[match.peptide] += 1
    aliases = _protein_aliases(protein_rows)
    return {
        name: _level_match(
            level,
            prepared[name],
            protein_rows,
            occurrences,
            site_counts,
            aliases,
            protein_group_separator,
        )
        for name, level in levels.items()
    }


def _level_match(
    level: PeptideLevelInput,
    peptides: tuple[str | None, ...],
    protein_rows: list[dict[str, object]],
    occurrences: Mapping[str, list[int]],
    site_counts: Mapping[str, int],
    aliases: Mapping[str, tuple[int, ...]],
    separator: str,
) -> PeptideLevelMatch:
    assignments = (
        level.frame.get_column(level.accession_column).to_list()
        if level.accession_column is not None
        else [None] * level.frame.height
    )
    rows: list[dict[str, object]] = []
    for peptide, assignment in zip(peptides, assignments, strict=True):
        record_indices = tuple(dict.fromkeys(occurrences.get(peptide or "", ())))
        reported = _members(assignment, separator) if level.accession_column is not None else None
        reported_records = _matching_records(reported or (), aliases)
        rows.append(
            {
                "peptide_in_fasta": bool(record_indices),
                "fasta_match_site_count": site_counts.get(peptide or "", 0),
                "fasta_matching_protein_count": len(record_indices),
                "fasta_matching_protein_ids": ";".join(
                    _required_text(protein_rows[index], "id") for index in record_indices
                ),
                "reported_member_count": None if reported is None else len(reported),
                "reported_members_in_fasta_count": (
                    None
                    if reported is None
                    else sum(bool(aliases.get(member)) for member in reported)
                ),
                "peptide_in_reported_protein": (
                    None
                    if reported is None
                    else bool(set(record_indices).intersection(reported_records))
                ),
            }
        )
    summary = pl.DataFrame(
        rows,
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
    matched_count = summary.get_column("peptide_in_fasta").sum() or 0
    return PeptideLevelMatch(
        summary=summary,
        coverage=PeptideCoverage(
            feature_count=len(peptides),
            unique_sequence_count=len({value for value in peptides if value is not None}),
            matched_feature_count=int(matched_count),
            unmatched_feature_count=len(peptides) - int(matched_count),
            match_site_count=sum(site_counts.get(value or "", 0) for value in set(peptides)),
        ),
    )


def _protein_aliases(rows: list[dict[str, object]]) -> dict[str, tuple[int, ...]]:
    collected: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        for column in ("id", "accession"):
            value = row.get(column)
            if isinstance(value, str) and value:
                collected[value].append(index)
    return {name: tuple(indices) for name, indices in collected.items()}


def _matching_records(
    members: tuple[str, ...],
    aliases: Mapping[str, tuple[int, ...]],
) -> set[int]:
    return {index for member in members for index in aliases.get(member, ())}


def _members(value: object, separator: str) -> tuple[str, ...]:
    if not isinstance(value, str):
        return ()
    return tuple(token.strip() for token in value.split(separator) if token.strip())


def _normalize_sequence(value: object, il_equivalent: bool) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    return _equivalent_sequence(value.strip().upper(), il_equivalent)


def _equivalent_sequence(value: str, enabled: bool) -> str:
    return value.replace("I", "L") if enabled else value


def _required_text(row: Mapping[str, object], column: str) -> str:
    value = row.get(column)
    if not isinstance(value, str):
        raise ValueError(f"protein frame column {column!r} contains a non-text value")
    return value
