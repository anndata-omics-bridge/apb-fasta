"""Thin command-line composition for independent APB2 FASTA operations."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Literal

from apb2.result_facade import read_parsed_levels, write_parsed_levels
from cyclopts import App
from loguru import logger
from protein_fasta.frame import ProteinDatabase, ProteinFormat, refseq, uniprotkb

from apb_fasta.annotation import FastaAnnotationParser, FastaAnnotationResult
from apb_fasta.configuration import FastaAnnotationParameters

app = App(
    name="apb-fasta",
    help="Verify peptides and merge protein annotations using FASTA files",
    help_on_error=True,
)

_FORMATS: dict[str, ProteinFormat] = {
    "refseq": refseq,
    "uniprotkb": uniprotkb,
}


@app.command
def verify_peptides(
    source: Path,
    *fasta_paths: Path,
    output: Path,
    formats: tuple[str, ...] = ("uniprotkb", "refseq"),
    backend: Literal["auto", "ahocorapy", "ahocorasick_rs"] = "auto",
    il_equivalent: bool = False,
    protein_group_separator: str = ";",
) -> int:
    """Verify stripped peptide sequences in SOURCE against FASTA_PATHS."""
    try:
        parser = _parser_for(
            source,
            fasta_paths,
            output=output,
            formats=formats,
            parameters=FastaAnnotationParameters(
                protein_group_separator=protein_group_separator,
                matcher_backend=backend,
                il_equivalent=il_equivalent,
            ),
        )
        result = parser.verify_peptides()
        _write_result(result, output)
        _report_peptide_verification(result)
    except (OSError, ValueError) as error:
        logger.error(str(error))
        return 1
    logger.info("wrote peptide-verified APB2 result to {}", output)
    return 0


@app.command
def merge_annotations(
    source: Path,
    *fasta_paths: Path,
    output: Path,
    formats: tuple[str, ...] = ("uniprotkb", "refseq"),
    protein_group_separator: str = ";",
) -> int:
    """Merge FASTA annotations into the reported protein groups in SOURCE."""
    try:
        parser = _parser_for(
            source,
            fasta_paths,
            output=output,
            formats=formats,
            parameters=FastaAnnotationParameters(
                protein_group_separator=protein_group_separator,
            ),
        )
        result = parser.merge_annotations()
        _write_result(result, output)
        _report_protein_annotations(result)
    except (OSError, ValueError) as error:
        logger.error(str(error))
        return 1
    logger.info("wrote FASTA-annotated APB2 result to {}", output)
    return 0


@app.command
def run(
    source: Path,
    *fasta_paths: Path,
    output: Path,
    formats: tuple[str, ...] = ("uniprotkb", "refseq"),
    backend: Literal["auto", "ahocorapy", "ahocorasick_rs"] = "auto",
    il_equivalent: bool = False,
    protein_group_separator: str = ";",
) -> int:
    """Verify peptides and merge protein annotations in one in-memory run."""
    try:
        parser = _parser_for(
            source,
            fasta_paths,
            output=output,
            formats=formats,
            parameters=FastaAnnotationParameters(
                protein_group_separator=protein_group_separator,
                matcher_backend=backend,
                il_equivalent=il_equivalent,
            ),
        )
        result = parser.run()
        _write_result(result, output)
        _report_peptide_verification(result)
        _report_protein_annotations(result)
    except (OSError, ValueError) as error:
        logger.error(str(error))
        return 1
    logger.info("wrote peptide-verified and FASTA-annotated APB2 result to {}", output)
    return 0


def _parser_for(
    source: Path,
    fasta_paths: tuple[Path, ...],
    /,
    *,
    output: Path,
    formats: tuple[str, ...],
    parameters: FastaAnnotationParameters,
) -> FastaAnnotationParser:
    if not fasta_paths:
        raise ValueError("at least one FASTA path is required")
    if output == source:
        raise ValueError("output must differ from source")
    if output.exists():
        raise ValueError(f"output already exists: {output}")
    unknown = tuple(name for name in formats if name not in _FORMATS)
    if unknown:
        raise ValueError(f"unknown protein format(s): {unknown}")
    protein_database = ProteinDatabase(*(_FORMATS[name] for name in formats))
    proteins = protein_database.parse(fasta_paths)
    return FastaAnnotationParser(
        read_parsed_levels(source),
        proteins,
        parameters=parameters,
    )


def _write_result(result: FastaAnnotationResult, output: Path) -> None:
    write_parsed_levels(result.parsed, output)


def _report_peptide_verification(result: FastaAnnotationResult) -> None:
    for level, coverage in result.reports.peptide_levels.items():
        logger.info(
            "level={} peptides_in_fasta={}/{} unmatched={} unique_sequences={} match_sites={}",
            level,
            coverage.matched_feature_count,
            coverage.feature_count,
            coverage.unmatched_feature_count,
            coverage.unique_sequence_count,
            coverage.match_site_count,
        )


def _report_protein_annotations(result: FastaAnnotationResult) -> None:
    coverage = result.reports.protein_groups
    if coverage is None:
        return
    logger.info(
        "protein_groups={} members_in_fasta={}/{} unmatched={} ambiguous={}",
        coverage.group_count,
        coverage.matched_member_count,
        coverage.member_count,
        coverage.unmatched_member_count,
        coverage.ambiguous_member_count,
    )


def main() -> int:
    """Run the console application."""
    result = app()
    return int(result) if result is not None else 0


if __name__ == "__main__":
    sys.exit(main())
