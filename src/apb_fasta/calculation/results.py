"""Immutable calculation values returned to the APB2 integration boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import polars as pl


@dataclass(frozen=True, slots=True)
class PeptideCoverage:
    """Aggregate FASTA coverage for one APB2 feature level."""

    feature_count: int
    unique_sequence_count: int
    matched_feature_count: int
    """Targets, rows apb2 did not mark ``apb_Decoy``, whose peptide the FASTA contains."""
    unmatched_feature_count: int
    """Targets whose peptide the FASTA does not contain."""
    match_site_count: int
    decoy_feature_count: int
    """Rows apb2 marked ``apb_Decoy``; a decoy is neither matched nor unmatched."""
    il_only_matched_feature_count: int
    """Targets in ``matched_feature_count`` whose peptide matches only as another I/L spelling."""


@dataclass(frozen=True, slots=True)
class PeptideLevelMatch:
    """Feature-aligned peptide facts and their aggregate coverage."""

    summary: pl.DataFrame
    coverage: PeptideCoverage


@dataclass(frozen=True, slots=True)
class ProteinGroupCoverage:
    """Aggregate FASTA coverage for reported protein-group members."""

    group_count: int
    member_count: int
    matched_member_count: int
    unmatched_member_count: int
    ambiguous_member_count: int


@dataclass(frozen=True, slots=True)
class ProteinGroupMatch:
    """Lossless member expansion, group summaries, and source relation."""

    members: pl.DataFrame
    member_key_columns: tuple[str, ...]
    summary: pl.DataFrame
    relation: pl.DataFrame
    coverage: ProteinGroupCoverage


@dataclass(frozen=True, slots=True)
class FastaAnnotationReports:
    """Persistable aggregate evidence from one annotation application."""

    peptide_levels: Mapping[str, PeptideCoverage]
    protein_groups: ProteinGroupCoverage | None
