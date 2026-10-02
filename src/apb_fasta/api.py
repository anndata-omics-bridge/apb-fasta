"""Public in-memory FASTA annotation API."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import version

import polars as pl
from apb2.api import ParsedLevels
from prozor.matching.automaton import resolve_backend

from apb_fasta.calculation.matching import match_peptide_levels
from apb_fasta.calculation.peptide_properties import peptide_level_properties
from apb_fasta.calculation.protein_groups import match_protein_groups
from apb_fasta.calculation.results import FastaAnnotationReports
from apb_fasta.configuration import (
    DEFAULT_FASTA_ANNOTATION_PARAMETERS,
    FastaAnnotationParameters,
)
from apb_fasta.errors import FastaAnnotationError
from apb_fasta.integration import (
    apply_peptide_matches,
    apply_peptide_properties,
    apply_protein_group_match,
    peptide_inputs,
    protein_frame_metadata,
    protein_group_input,
    validate_protein_frame,
)

_NO_PEPTIDE_LEVEL = (
    "result contains no peptide-derived level with canonical ProForma_peptide values"
)


@dataclass(frozen=True, slots=True)
class FastaAnnotationResult:
    """A replacement APB2 result and its FASTA reports."""

    parsed: ParsedLevels
    reports: FastaAnnotationReports


class FastaAnnotator:
    """Bind one protein database and its FASTA annotation behavior."""

    __slots__ = ("_parameters", "_proteins")

    def __init__(
        self,
        proteins: pl.DataFrame,
        /,
        *,
        parameters: FastaAnnotationParameters = DEFAULT_FASTA_ANNOTATION_PARAMETERS,
    ) -> None:
        """Validate and bind a reusable protein database and configuration."""
        validate_protein_frame(proteins)
        self._proteins = proteins
        self._parameters = parameters

    def verify_peptides(self, parsed: ParsedLevels, /) -> FastaAnnotationResult:
        """Verify every canonical stripped peptide against the protein sequences."""
        inputs = peptide_inputs(parsed)
        if not inputs:
            raise FastaAnnotationError(_NO_PEPTIDE_LEVEL)
        peptide_levels = match_peptide_levels(
            inputs,
            self._proteins,
            backend=self._parameters.matcher_backend,
            il_equivalent=self._parameters.il_equivalent,
            protein_group_separator=self._parameters.protein_group_separator,
        )
        replacement, reports = apply_peptide_matches(
            parsed,
            peptide_levels,
            requested_backend=self._parameters.matcher_backend,
            resolved_backend=resolve_backend(self._parameters.matcher_backend),
            il_equivalent=self._parameters.il_equivalent,
            protein_metadata=protein_frame_metadata(self._proteins),
        )
        return FastaAnnotationResult(parsed=replacement, reports=reports)

    def merge_annotations(self, parsed: ParsedLevels, /) -> FastaAnnotationResult:
        """Merge FASTA annotations for every reported protein-group member."""
        protein_input = protein_group_input(parsed)
        if protein_input is None:
            raise FastaAnnotationError("result contains no protein level to annotate")
        protein_groups = match_protein_groups(
            protein_input,
            self._proteins,
            separator=self._parameters.protein_group_separator,
        )
        replacement, reports = apply_protein_group_match(
            parsed,
            protein_groups,
            protein_group_separator=self._parameters.protein_group_separator,
            protein_metadata=protein_frame_metadata(self._proteins),
        )
        return FastaAnnotationResult(parsed=replacement, reports=reports)

    def annotate(self, parsed: ParsedLevels, /) -> FastaAnnotationResult:
        """Verify peptides and then merge protein annotations in memory."""
        verified = self.verify_peptides(parsed)
        annotated = self.merge_annotations(verified.parsed)
        return FastaAnnotationResult(
            parsed=annotated.parsed,
            reports=FastaAnnotationReports(
                peptide_levels=verified.reports.peptide_levels,
                protein_groups=annotated.reports.protein_groups,
            ),
        )


def add_peptide_properties(parsed: ParsedLevels, /) -> ParsedLevels:
    """Attach sequence-derived properties to every peptide-derived level.

    Writes a feature-aligned ``varm["peptide_properties"]`` computed by
    ``protein_fasta.peptide_frame.peptide_property_frame`` from each level's ``ProForma_peptide``.
    No protein database is needed. A feature without a sequence, or with a residue outside the 20
    standard amino acids, has null properties.

    Raises:
        FastaAnnotationError: If the result has no peptide-derived level, already contains the
            output, or already records the operation.
        ValueError: If a sequence is not stripped upper-case letters.
    """
    inputs = peptide_inputs(parsed)
    if not inputs:
        raise FastaAnnotationError(_NO_PEPTIDE_LEVEL)
    properties = {
        name: peptide_level_properties(level.frame, level.sequence_column)
        for name, level in inputs.items()
    }
    return apply_peptide_properties(
        parsed, properties, protein_fasta_version=version("protein-fasta")
    )


__all__ = [
    "FastaAnnotationResult",
    "FastaAnnotator",
    "add_peptide_properties",
]
