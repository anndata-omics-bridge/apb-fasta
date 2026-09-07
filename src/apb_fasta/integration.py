"""Translate between APB2 result values and FASTA calculations."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from typing import cast

import polars as pl
from apb2.result_facade import (
    AnnotationTable,
    FeatureRelation,
    JsonValue,
    ParsedLevelName,
    ParsedLevels,
)

from apb_fasta.calculation.matching import PeptideLevelInput
from apb_fasta.calculation.protein_groups import ProteinGroupInput
from apb_fasta.calculation.results import (
    FastaAnnotationReports,
    PeptideLevelMatch,
    ProteinGroupMatch,
)
from apb_fasta.errors import FastaAnnotationError

FASTA_VALIDATION_NAME = "fasta_validation"
FASTA_SUMMARY_NAME = "fasta"
MEMBER_TABLE_NAME = "fasta_protein_group_members"
MEMBER_RELATION_NAME = "fasta_protein_group_membership"
PEPTIDE_VERIFICATION_OPERATION = "peptide_verification"
PROTEIN_ANNOTATION_OPERATION = "protein_annotation"
_PEPTIDE_COLUMN = "ProForma_peptide"
_RESERVED_MEMBER_COLUMNS = frozenset({"source_level", "member_ordinal", "match_ordinal"})


def validate_protein_frame(proteins: pl.DataFrame, /) -> None:
    """Require the stable protein-fasta columns used by both calculations."""
    required = {
        "id",
        "description",
        "sequence",
        "fasta_source_path",
        "fasta_source_checksum",
        "fasta_source_ordinal",
        "fasta_record_ordinal",
    }
    missing = sorted(required.difference(proteins.columns))
    if missing:
        raise FastaAnnotationError(f"protein frame is missing required columns {missing}")
    if proteins.get_column("id").null_count():
        raise FastaAnnotationError("protein frame contains a missing id")
    if proteins.get_column("sequence").null_count():
        raise FastaAnnotationError("protein frame contains a missing sequence")


def peptide_inputs(parsed: ParsedLevels, /) -> dict[str, PeptideLevelInput]:
    """Extract every peptide-derived level carrying APB2's canonical sequence."""
    result: dict[str, PeptideLevelInput] = {}
    for name, level in parsed.levels.items():
        if name == "protein" or _PEPTIDE_COLUMN not in level.var.frame.columns:
            continue
        result[name] = PeptideLevelInput(
            frame=level.var.frame,
            sequence_column=_PEPTIDE_COLUMN,
            accession_column=_fasta_accession_column(level.uns, level.var.frame),
        )
    return result


def protein_group_input(parsed: ParsedLevels, /) -> ProteinGroupInput | None:
    """Extract the protein axis and its persisted FASTA-accession role."""
    level = parsed.levels.get("protein")
    if level is None:
        return None
    accession_column = _fasta_accession_column(level.uns, level.var.frame)
    if accession_column is None:
        raise FastaAnnotationError("protein level does not declare column_roles.fasta_accessions")
    collisions = _RESERVED_MEMBER_COLUMNS.intersection(level.var.key_columns)
    if collisions:
        raise FastaAnnotationError(
            f"protein variable keys collide with FASTA member columns {sorted(collisions)}"
        )
    return ProteinGroupInput(
        frame=level.var.frame,
        key_columns=level.var.key_columns,
        accession_column=accession_column,
    )


def validate_peptide_output_names(
    parsed: ParsedLevels,
    matches: dict[str, PeptideLevelMatch],
    /,
) -> None:
    """Reject peptide-verification collisions before constructing a replacement."""
    _validate_operation_metadata(parsed, PEPTIDE_VERIFICATION_OPERATION)
    for name in matches:
        level_name = cast(ParsedLevelName, name)
        if FASTA_VALIDATION_NAME in parsed.levels[level_name].varm:
            raise FastaAnnotationError(
                f"level {name!r} already contains varm[{FASTA_VALIDATION_NAME!r}]"
            )
        if "fasta" in parsed.levels[level_name].metadata:
            raise FastaAnnotationError(f"level {name!r} metadata already contains 'fasta'")


def validate_protein_output_names(parsed: ParsedLevels, /) -> None:
    """Reject protein-annotation collisions before constructing a replacement."""
    _validate_operation_metadata(parsed, PROTEIN_ANNOTATION_OPERATION)
    protein = parsed.levels["protein"]
    if FASTA_SUMMARY_NAME in protein.varm:
        raise FastaAnnotationError(f"level 'protein' already contains varm[{FASTA_SUMMARY_NAME!r}]")
    if MEMBER_TABLE_NAME in parsed.annotation_tables:
        raise FastaAnnotationError(
            f"result already contains annotation table {MEMBER_TABLE_NAME!r}"
        )
    if MEMBER_RELATION_NAME in parsed.feature_relations:
        raise FastaAnnotationError(
            f"result already contains feature relation {MEMBER_RELATION_NAME!r}"
        )


def apply_peptide_matches(
    parsed: ParsedLevels,
    matches: dict[str, PeptideLevelMatch],
    /,
    *,
    requested_backend: str,
    resolved_backend: str,
    il_equivalent: bool,
    protein_metadata: dict[str, JsonValue],
) -> tuple[ParsedLevels, FastaAnnotationReports]:
    """Attach peptide verification to a deep replacement of the APB2 result."""
    validate_peptide_output_names(parsed, matches)
    result = deepcopy(parsed)
    for name, match in matches.items():
        level_name = cast(ParsedLevelName, name)
        result.levels[level_name].varm[FASTA_VALIDATION_NAME] = match.summary.clone()
        level_metadata = result.levels[level_name].metadata
        level_metadata["fasta"] = {
            PEPTIDE_VERIFICATION_OPERATION: {
                **cast(dict[str, JsonValue], asdict(match.coverage)),
                "requested_backend": requested_backend,
                "resolved_backend": resolved_backend,
                **protein_metadata,
            }
        }
    _record_operation_metadata(
        result,
        PEPTIDE_VERIFICATION_OPERATION,
        {
            "requested_backend": requested_backend,
            "resolved_backend": resolved_backend,
            "il_equivalent": il_equivalent,
            **protein_metadata,
        },
    )
    reports = FastaAnnotationReports(
        peptide_levels={name: match.coverage for name, match in matches.items()},
        protein_groups=None,
    )
    return result, reports


def apply_protein_group_match(
    parsed: ParsedLevels,
    match: ProteinGroupMatch,
    /,
    *,
    protein_group_separator: str,
    protein_metadata: dict[str, JsonValue],
) -> tuple[ParsedLevels, FastaAnnotationReports]:
    """Attach protein FASTA annotations to a deep replacement of the APB2 result."""
    validate_protein_output_names(parsed)
    result = deepcopy(parsed)
    result.levels["protein"].varm[FASTA_SUMMARY_NAME] = match.summary.clone()
    result.annotation_tables[MEMBER_TABLE_NAME] = AnnotationTable(
        frame=match.members.clone(),
        key_columns=match.member_key_columns,
        metadata={"producer": "apb-fasta", "schema_version": "1"},
    )
    result.feature_relations[MEMBER_RELATION_NAME] = FeatureRelation(
        annotation_table=MEMBER_TABLE_NAME,
        target_level="protein",
        coordinates=match.relation.clone(),
        metadata={"producer": "apb-fasta", "semantic": "member_of"},
    )
    _record_operation_metadata(
        result,
        PROTEIN_ANNOTATION_OPERATION,
        {
            "protein_group_separator": protein_group_separator,
            **protein_metadata,
        },
    )
    return result, FastaAnnotationReports(
        peptide_levels={},
        protein_groups=match.coverage,
    )


def protein_frame_metadata(proteins: pl.DataFrame, /) -> dict[str, JsonValue]:
    """Return bounded source and schema provenance for one protein frame."""
    sources = (
        proteins.select(
            "fasta_source_ordinal",
            "fasta_source_path",
            "fasta_source_checksum",
        )
        .unique(maintain_order=True)
        .sort("fasta_source_ordinal")
        .to_dicts()
    )
    return {
        "protein_count": proteins.height,
        "protein_columns": list(proteins.columns),
        "sources": {
            str(_json_scalar(source["fasta_source_ordinal"])): {
                "ordinal": _json_scalar(source["fasta_source_ordinal"]),
                "path": _json_scalar(source["fasta_source_path"]),
                "checksum": _json_scalar(source["fasta_source_checksum"]),
            }
            for source in sources
        },
    }


def _validate_operation_metadata(parsed: ParsedLevels, operation: str) -> None:
    metadata = parsed.metadata.get("fasta")
    if metadata is None:
        return
    if not isinstance(metadata, dict):
        raise FastaAnnotationError("result metadata 'fasta' section is not an object")
    if operation in metadata:
        raise FastaAnnotationError(f"result already contains FASTA operation {operation!r}")


def _record_operation_metadata(
    parsed: ParsedLevels,
    operation: str,
    operation_metadata: dict[str, JsonValue],
) -> None:
    existing = parsed.metadata.get("fasta")
    metadata: dict[str, JsonValue]
    if existing is None:
        metadata = {"schema_version": "1"}
    else:
        metadata = cast(dict[str, JsonValue], deepcopy(existing))
    metadata[operation] = operation_metadata
    parsed.metadata["fasta"] = metadata


def _fasta_accession_column(
    uns: dict[str, JsonValue],
    frame: pl.DataFrame,
) -> str | None:
    roles = uns.get("column_roles")
    if not isinstance(roles, dict):
        return None
    value = roles.get("fasta_accessions")
    return value if isinstance(value, str) and value in frame.columns else None


def _json_scalar(value: object) -> bool | int | float | str | None:
    if value is None or isinstance(value, bool | int | float | str):
        return value
    raise FastaAnnotationError(f"protein source provenance contains {type(value).__name__}")
