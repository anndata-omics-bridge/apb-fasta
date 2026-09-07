"""Dataset- and protein-bound FASTA operations."""

from __future__ import annotations

from dataclasses import dataclass, field

import polars as pl
from apb2.result_facade import ParsedLevels
from prozor.matching.automaton import resolve_backend

from apb_fasta.calculation.matching import match_peptide_levels
from apb_fasta.calculation.protein_groups import match_protein_groups
from apb_fasta.calculation.results import FastaAnnotationReports
from apb_fasta.configuration import (
    DEFAULT_FASTA_ANNOTATION_PARAMETERS,
    FastaAnnotationParameters,
)
from apb_fasta.errors import FastaAnnotationError
from apb_fasta.integration import (
    apply_peptide_matches,
    apply_protein_group_match,
    peptide_inputs,
    protein_frame_metadata,
    protein_group_input,
    validate_protein_frame,
)


@dataclass(frozen=True, slots=True)
class FastaAnnotationResult:
    """A replacement APB2 result and its FASTA reports."""

    parsed: ParsedLevels
    reports: FastaAnnotationReports


@dataclass(frozen=True, slots=True)
class FastaAnnotationParser:
    """Apply independent FASTA operations to one APB2 result and protein frame."""

    anndata: ParsedLevels
    proteins: pl.DataFrame
    parameters: FastaAnnotationParameters = field(
        default=DEFAULT_FASTA_ANNOTATION_PARAMETERS,
        kw_only=True,
    )

    def verify_peptides(self) -> FastaAnnotationResult:
        """Verify every canonical stripped peptide against the protein sequences."""
        validate_protein_frame(self.proteins)
        inputs = peptide_inputs(self.anndata)
        if not inputs:
            raise FastaAnnotationError(
                "result contains no peptide-derived level with canonical ProForma_peptide values"
            )
        peptide_levels = match_peptide_levels(
            inputs,
            self.proteins,
            backend=self.parameters.matcher_backend,
            il_equivalent=self.parameters.il_equivalent,
            protein_group_separator=self.parameters.protein_group_separator,
        )
        parsed, reports = apply_peptide_matches(
            self.anndata,
            peptide_levels,
            requested_backend=self.parameters.matcher_backend,
            resolved_backend=resolve_backend(self.parameters.matcher_backend),
            il_equivalent=self.parameters.il_equivalent,
            protein_metadata=protein_frame_metadata(self.proteins),
        )
        return FastaAnnotationResult(parsed=parsed, reports=reports)

    def merge_annotations(self) -> FastaAnnotationResult:
        """Merge FASTA annotations for every reported protein-group member."""
        validate_protein_frame(self.proteins)
        protein_input = protein_group_input(self.anndata)
        if protein_input is None:
            raise FastaAnnotationError("result contains no protein level to annotate")
        protein_groups = match_protein_groups(
            protein_input,
            self.proteins,
            separator=self.parameters.protein_group_separator,
        )
        parsed, reports = apply_protein_group_match(
            self.anndata,
            protein_groups,
            protein_group_separator=self.parameters.protein_group_separator,
            protein_metadata=protein_frame_metadata(self.proteins),
        )
        return FastaAnnotationResult(parsed=parsed, reports=reports)

    def run(self) -> FastaAnnotationResult:
        """Verify peptides and then merge protein annotations in memory."""
        verified = self.verify_peptides()
        annotated = FastaAnnotationParser(
            verified.parsed,
            self.proteins,
            parameters=self.parameters,
        ).merge_annotations()
        return FastaAnnotationResult(
            parsed=annotated.parsed,
            reports=FastaAnnotationReports(
                peptide_levels=verified.reports.peptide_levels,
                protein_groups=annotated.reports.protein_groups,
            ),
        )
