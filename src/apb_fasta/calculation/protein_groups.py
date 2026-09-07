"""Protein-group expansion and FASTA annotation over Polars values."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass

import polars as pl

from apb_fasta.calculation.results import ProteinGroupCoverage, ProteinGroupMatch

_OPTIONAL_FASTA_COLUMNS = (
    "is_decoy",
    "is_contaminant",
    "accession",
    "protein_name",
    "gene_name",
    "organism_name",
    "taxonomy_id",
    "database",
    "review_status",
    "fasta_source_path",
    "fasta_source_checksum",
    "fasta_source_ordinal",
    "fasta_record_ordinal",
)


@dataclass(frozen=True, slots=True)
class ProteinGroupInput:
    """The protein feature axis and its declared FASTA accession column."""

    frame: pl.DataFrame
    key_columns: tuple[str, ...]
    accession_column: str


def match_protein_groups(
    source: ProteinGroupInput,
    proteins: pl.DataFrame,
    /,
    *,
    separator: str,
) -> ProteinGroupMatch:
    """Expand every reported member and retain every matching FASTA record."""
    protein_rows = proteins.to_dicts()
    aliases = _protein_aliases(protein_rows)
    members_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    relation_rows: list[dict[str, object]] = []
    matched_members = 0
    unmatched_members = 0
    ambiguous_members = 0
    member_count = 0
    for group_index, group in enumerate(source.frame.to_dicts()):
        raw_group = group.get(source.accession_column)
        members = _members(raw_group, separator)
        group_matched = 0
        group_unmatched = 0
        group_ambiguous = 0
        descriptions: list[str] = []
        gene_names: list[str] = []
        organism_names: list[str] = []
        for member_ordinal, member in enumerate(members):
            member_count += 1
            matches = aliases.get(member, ())
            if not matches:
                unmatched_members += 1
                group_unmatched += 1
                members_rows.append(
                    _member_row(
                        source,
                        group,
                        raw_group,
                        member,
                        member_ordinal,
                        0,
                        "unmatched",
                        None,
                    )
                )
                relation_rows.append(
                    {"row": len(members_rows) - 1, "column": group_index, "value": 1.0}
                )
                continue
            matched_members += 1
            group_matched += 1
            status = "ambiguous" if len(matches) > 1 else "matched"
            if len(matches) > 1:
                ambiguous_members += 1
                group_ambiguous += 1
            for match_ordinal, protein_index in enumerate(matches):
                protein = protein_rows[protein_index]
                members_rows.append(
                    _member_row(
                        source,
                        group,
                        raw_group,
                        member,
                        member_ordinal,
                        match_ordinal,
                        status,
                        protein,
                    )
                )
                relation_rows.append(
                    {"row": len(members_rows) - 1, "column": group_index, "value": 1.0}
                )
                _append_text(descriptions, protein.get("description"))
                _append_text(gene_names, protein.get("gene_name"))
                _append_text(organism_names, protein.get("organism_name"))
        summaries.append(
            {
                "reported_member_count": len(members),
                "matched_member_count": group_matched,
                "unmatched_member_count": group_unmatched,
                "ambiguous_member_count": group_ambiguous,
                "all_members_in_fasta": bool(members) and group_unmatched == 0,
                "any_member_in_fasta": group_matched > 0,
                "fasta_descriptions": ";".join(descriptions),
                "fasta_gene_names": ";".join(gene_names),
                "fasta_organism_names": ";".join(organism_names),
            }
        )
    member_keys = ("source_level", *source.key_columns, "member_ordinal", "match_ordinal")
    return ProteinGroupMatch(
        members=pl.DataFrame(members_rows, schema=_member_schema(source, proteins)),
        member_key_columns=member_keys,
        summary=pl.DataFrame(summaries, schema=_summary_schema()),
        relation=pl.DataFrame(
            relation_rows,
            schema={"row": pl.Int64, "column": pl.Int64, "value": pl.Float64},
        ),
        coverage=ProteinGroupCoverage(
            group_count=source.frame.height,
            member_count=member_count,
            matched_member_count=matched_members,
            unmatched_member_count=unmatched_members,
            ambiguous_member_count=ambiguous_members,
        ),
    )


def _member_row(
    source: ProteinGroupInput,
    group: Mapping[str, object],
    raw_group: object,
    member: str,
    member_ordinal: int,
    match_ordinal: int,
    status: str,
    protein: Mapping[str, object] | None,
) -> dict[str, object]:
    sequence = None if protein is None else protein.get("sequence")
    row: dict[str, object] = {
        "source_level": "protein",
        **{name: group[name] for name in source.key_columns},
        "member_ordinal": member_ordinal,
        "match_ordinal": match_ordinal,
        "protein_group": raw_group,
        "protein_member": member,
        "normalized_matching_key": member,
        "match_status": status,
        "fasta_id": None if protein is None else protein.get("id"),
        "fasta_description": None if protein is None else protein.get("description"),
        "fasta_sequence_length": (None if not isinstance(sequence, str) else len(sequence)),
    }
    for name in _OPTIONAL_FASTA_COLUMNS:
        row[name] = None if protein is None else protein.get(name)
    return row


def _member_schema(
    source: ProteinGroupInput,
    proteins: pl.DataFrame,
) -> dict[str, pl.DataType | type[pl.DataType]]:
    schema: dict[str, pl.DataType | type[pl.DataType]] = {
        "source_level": pl.String,
        **{name: source.frame.schema[name] for name in source.key_columns},
        "member_ordinal": pl.Int64,
        "match_ordinal": pl.Int64,
        "protein_group": source.frame.schema[source.accession_column],
        "protein_member": pl.String,
        "normalized_matching_key": pl.String,
        "match_status": pl.String,
        "fasta_id": pl.String,
        "fasta_description": pl.String,
        "fasta_sequence_length": pl.UInt64,
    }
    for name in _OPTIONAL_FASTA_COLUMNS:
        schema[name] = proteins.schema.get(name, pl.Null)
    return schema


def _summary_schema() -> dict[str, type[pl.DataType]]:
    return {
        "reported_member_count": pl.UInt64,
        "matched_member_count": pl.UInt64,
        "unmatched_member_count": pl.UInt64,
        "ambiguous_member_count": pl.UInt64,
        "all_members_in_fasta": pl.Boolean,
        "any_member_in_fasta": pl.Boolean,
        "fasta_descriptions": pl.String,
        "fasta_gene_names": pl.String,
        "fasta_organism_names": pl.String,
    }


def _protein_aliases(rows: list[dict[str, object]]) -> dict[str, tuple[int, ...]]:
    collected: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        aliases = {
            value
            for name in ("id", "accession")
            if isinstance((value := row.get(name)), str) and value
        }
        for alias in aliases:
            collected[alias].append(index)
    return {name: tuple(indices) for name, indices in collected.items()}


def _members(value: object, separator: str) -> tuple[str, ...]:
    if not isinstance(value, str):
        return ()
    return tuple(token.strip() for token in value.split(separator) if token.strip())


def _append_text(values: list[str], value: object) -> None:
    if isinstance(value, str) and value and value not in values:
        values.append(value)
