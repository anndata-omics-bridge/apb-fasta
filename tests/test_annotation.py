"""Public FASTA annotation lifecycle and APB2 persistence."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from apb2.api import ParsedLevel, ParsedLevels, read_parsed_levels, write_parsed_levels

from apb_fasta.api import FastaAnnotator, add_peptide_properties
from apb_fasta.cli import app
from apb_fasta.errors import FastaAnnotationError


def _level(
    name: str,
    var: pl.DataFrame,
    key_columns: tuple[str, ...],
    column_roles: dict[str, str],
) -> ParsedLevel:
    return ParsedLevel.build(
        pl.DataFrame({"Run": ["run1"]}),
        ("Run",),
        var,
        key_columns,
        column_roles,
        primary_layer="Intensity",
        abundance={
            "Intensity": pl.DataFrame({"obs_0": [float(index + 1) for index in range(var.height)]})
        },
        uns={"quantification_level": name},
    )


def _parsed() -> ParsedLevels:
    peptide = _level(
        "peptide",
        pl.DataFrame(
            {
                "ProForma_peptide": ["PEPTIDE", "OTHER", "MISSING"],
                "Proteins": ["P1;P2", "P2", "P9"],
            }
        ),
        ("ProForma_peptide",),
        {"fasta_accessions": "Proteins"},
    )
    protein = _level(
        "protein",
        pl.DataFrame(
            {
                "ProteinGroup": ["P1;P2", "P3", "P9"],
                "Proteins": ["P1;P2", "P3", "P9"],
            }
        ),
        ("ProteinGroup",),
        {"fasta_accessions": "Proteins"},
    )
    return ParsedLevels(levels={"peptide": peptide, "protein": protein}, uns={})


def _proteins() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "id": ["sp|P1|ONE", "sp|P2|TWO", "REV_sp|P2|TWO"],
            "description": ["Protein one", "Protein two", "Decoy protein two"],
            "sequence": ["MPEPTIDEK", "XXPEPTIDEXX", "MOTHERK"],
            "accession": ["P1", "P2", "P2"],
            "protein_name": ["One", "Two", "Two decoy"],
            "gene_name": ["G1", "G2", "G2"],
            "organism_name": ["Human", "Human", "Human"],
            "is_decoy": [False, False, True],
            "is_contaminant": [False, False, False],
            "fasta_source_path": ["human.fasta", "human.fasta", "decoy.fasta"],
            "fasta_source_checksum": ["a", "a", "b"],
            "fasta_source_ordinal": [0, 0, 1],
            "fasta_record_ordinal": [0, 1, 0],
        }
    )


def test_annotate_verifies_peptides_and_merges_annotations_without_mutating_input() -> None:
    source = _parsed()

    result = FastaAnnotator(_proteins()).annotate(source)

    assert source.levels["peptide"].varm == {}
    assert source.levels["protein"].varm == {}
    assert "fasta_validation" in result.parsed.levels["peptide"].varm
    assert "fasta" in result.parsed.levels["protein"].varm
    assert "fasta_protein_group_members" in result.parsed.annotation_tables
    assert "fasta_protein_group_membership" in result.parsed.feature_relations
    assert result.reports.peptide_levels["peptide"].matched_feature_count == 2
    assert result.reports.protein_groups is not None
    assert result.reports.protein_groups.member_count == 4
    assert result.reports.protein_groups.matched_member_count == 2
    assert result.reports.protein_groups.ambiguous_member_count == 1
    metadata = result.parsed.metadata["fasta"]
    assert isinstance(metadata, dict)
    assert set(metadata) == {"schema_version", "provenance"}
    assert metadata["schema_version"] == "2"
    assert isinstance(metadata["provenance"], dict)
    assert set(metadata["provenance"]) == {"peptide_verification", "protein_annotation"}


def test_operations_run_independently() -> None:
    parsed = _parsed()
    annotator = FastaAnnotator(_proteins())

    verified = annotator.verify_peptides(parsed)
    annotated = annotator.merge_annotations(parsed)

    assert "fasta_validation" in verified.parsed.levels["peptide"].varm
    assert verified.parsed.levels["protein"].varm == {}
    assert verified.parsed.annotation_tables == {}
    assert verified.reports.protein_groups is None
    assert annotated.parsed.levels["peptide"].varm == {}
    assert "fasta" in annotated.parsed.levels["protein"].varm
    assert annotated.reports.peptide_levels == {}
    assert annotated.reports.protein_groups is not None


def test_operations_compose_in_memory() -> None:
    annotator = FastaAnnotator(_proteins())
    verified = annotator.verify_peptides(_parsed())

    complete = annotator.merge_annotations(verified.parsed)

    assert "fasta_validation" in complete.parsed.levels["peptide"].varm
    assert "fasta" in complete.parsed.levels["protein"].varm
    metadata = complete.parsed.metadata["fasta"]
    assert isinstance(metadata, dict)
    assert set(metadata) == {"schema_version", "provenance"}
    assert isinstance(metadata["provenance"], dict)
    assert set(metadata["provenance"]) == {"peptide_verification", "protein_annotation"}


def test_annotate_equals_explicit_composition() -> None:
    parsed = _parsed()
    annotator = FastaAnnotator(_proteins())

    combined = annotator.annotate(parsed)
    verified = annotator.verify_peptides(parsed)
    explicit = annotator.merge_annotations(verified.parsed)

    assert combined.reports.peptide_levels == verified.reports.peptide_levels
    assert combined.reports.protein_groups == explicit.reports.protein_groups
    assert combined.parsed.metadata == explicit.parsed.metadata
    assert (
        combined.parsed.levels["peptide"]
        .varm["fasta_validation"]
        .equals(explicit.parsed.levels["peptide"].varm["fasta_validation"])
    )
    assert (
        combined.parsed.levels["protein"]
        .varm["fasta"]
        .equals(explicit.parsed.levels["protein"].varm["fasta"])
    )
    combined_members = combined.parsed.annotation_tables["fasta_protein_group_members"]
    explicit_members = explicit.parsed.annotation_tables["fasta_protein_group_members"]
    assert combined_members.frame.equals(explicit_members.frame)
    combined_relation = combined.parsed.feature_relations["fasta_protein_group_membership"]
    explicit_relation = explicit.parsed.feature_relations["fasta_protein_group_membership"]
    assert combined_relation.coordinates.equals(explicit_relation.coordinates)


def test_one_annotator_processes_independent_results() -> None:
    annotator = FastaAnnotator(_proteins())
    first_input = _parsed()
    second_input = _parsed()

    first = annotator.verify_peptides(first_input)
    second = annotator.verify_peptides(second_input)

    assert first.parsed is not second.parsed
    assert (
        first.parsed.levels["peptide"]
        .varm["fasta_validation"]
        .equals(second.parsed.levels["peptide"].varm["fasta_validation"])
    )
    assert first_input.levels["peptide"].varm == {}
    assert second_input.levels["peptide"].varm == {}


def test_operations_compose_in_reverse_order() -> None:
    annotator = FastaAnnotator(_proteins())
    annotated = annotator.merge_annotations(_parsed())

    complete = annotator.verify_peptides(annotated.parsed)

    assert "fasta_validation" in complete.parsed.levels["peptide"].varm
    assert "fasta" in complete.parsed.levels["protein"].varm


@pytest.mark.parametrize("suffix", [".h5mu", ".parquet", ".duckdb"])
def test_annotated_result_round_trips_through_multilevel_formats(
    suffix: str,
    tmp_path: Path,
) -> None:
    result = FastaAnnotator(_proteins()).annotate(_parsed()).parsed
    target = tmp_path / f"annotated{suffix}"

    write_parsed_levels(result, target)
    restored = read_parsed_levels(target)

    assert restored.levels["peptide"].varm["fasta_validation"].height == 3
    members = restored.annotation_tables["fasta_protein_group_members"]
    assert members.frame.height == 5
    relation = restored.feature_relations["fasta_protein_group_membership"]
    assert relation.coordinates["column"].to_list() == [0, 0, 0, 1, 2]


def test_peptide_only_annotation_round_trips_through_h5ad(tmp_path: Path) -> None:
    parsed = _parsed()
    peptide_only = ParsedLevels(levels={"peptide": parsed.levels["peptide"]}, uns={})
    result = FastaAnnotator(_proteins()).verify_peptides(peptide_only).parsed
    target = tmp_path / "annotated.h5ad"

    write_parsed_levels(result, target)
    restored = read_parsed_levels(target)

    assert restored.levels["peptide"].varm["fasta_validation"].height == 3
    assert restored.annotation_tables == {}
    assert restored.feature_relations == {}


def test_annotator_rejects_owned_output_collisions() -> None:
    parsed = _parsed()
    parsed.levels["peptide"].varm["fasta_validation"] = pl.DataFrame({"old": [1, 2, 3]})

    with pytest.raises(FastaAnnotationError, match="already contains"):
        FastaAnnotator(_proteins()).verify_peptides(parsed)


def test_annotator_rejects_repeating_an_applied_operation() -> None:
    annotator = FastaAnnotator(_proteins())
    verified = annotator.verify_peptides(_parsed())

    with pytest.raises(FastaAnnotationError, match="peptide_verification"):
        annotator.verify_peptides(verified.parsed)


def test_constructor_requires_protein_fasta_database_schema() -> None:
    with pytest.raises(FastaAnnotationError, match="missing required columns"):
        FastaAnnotator(_proteins().drop("fasta_source_checksum"))


@pytest.mark.parametrize(
    ("command", "has_verification", "has_annotations"),
    [
        ("verify-peptides", True, False),
        ("merge-annotations", False, True),
        ("run", True, True),
    ],
)
def test_cli_exposes_independent_and_combined_operations(
    command: str,
    has_verification: bool,
    has_annotations: bool,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.h5mu"
    fasta = tmp_path / "proteins.fasta"
    output = tmp_path / "annotated.h5mu"
    write_parsed_levels(_parsed(), source)
    fasta.write_text(
        (
            ">sp|P1|ONE Protein one OS=Human OX=9606 PE=1 SV=1\nMPEPTIDEK\n"
            ">sp|P2|TWO Protein two OS=Human OX=9606 PE=1 SV=1\nMOTHERK\n"
        ),
        encoding="utf-8",
    )

    status = app(
        [command, str(source), str(fasta), "--output", str(output)],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    restored = read_parsed_levels(output)
    assert ("fasta_validation" in restored.levels["peptide"].varm) is has_verification
    assert ("fasta" in restored.levels["protein"].varm) is has_annotations


def test_peptide_properties_attach_to_peptide_levels_without_mutating_input() -> None:
    source = _parsed()

    result = add_peptide_properties(source)

    assert source.levels["peptide"].varm == {}
    assert "fasta" not in source.metadata
    assert result.levels["protein"].varm == {}
    properties = result.levels["peptide"].varm["peptide_properties"]
    # "OTHER" contains pyrrolysine (O), so its properties are null.
    assert properties.get_column("length").to_list() == [7, None, 7]
    metadata = result.metadata["fasta"]
    assert isinstance(metadata, dict)
    provenance = metadata["provenance"]
    assert isinstance(provenance, dict)
    operation = provenance["peptide_properties"]
    assert isinstance(operation, dict)
    assert isinstance(operation["protein_fasta_version"], str)


def test_peptide_properties_compose_with_fasta_annotation() -> None:
    annotated = FastaAnnotator(_proteins()).annotate(_parsed()).parsed

    complete = add_peptide_properties(annotated)

    assert set(complete.levels["peptide"].varm) == {"fasta_validation", "peptide_properties"}
    metadata = complete.metadata["fasta"]
    assert isinstance(metadata, dict)
    assert isinstance(metadata["provenance"], dict)
    assert set(metadata["provenance"]) == {
        "peptide_verification",
        "protein_annotation",
        "peptide_properties",
    }


def test_peptide_properties_refuse_a_second_application() -> None:
    once = add_peptide_properties(_parsed())

    with pytest.raises(FastaAnnotationError, match="already contains"):
        add_peptide_properties(once)


def test_peptide_properties_require_a_peptide_level() -> None:
    protein_only = ParsedLevels(levels={"protein": _parsed().levels["protein"]}, uns={})

    with pytest.raises(FastaAnnotationError, match="no peptide-derived level"):
        add_peptide_properties(protein_only)


@pytest.mark.parametrize("suffix", [".h5mu", ".parquet", ".duckdb"])
def test_peptide_properties_round_trip_through_multilevel_formats(
    suffix: str,
    tmp_path: Path,
) -> None:
    result = add_peptide_properties(_parsed())
    target = tmp_path / f"properties{suffix}"

    write_parsed_levels(result, target)
    restored = read_parsed_levels(target)

    properties = restored.levels["peptide"].varm["peptide_properties"]
    assert properties.get_column("length").to_list() == [7, None, 7]
    assert properties.get_column("c_terminal_residue").to_list() == ["E", None, "G"]
    assert properties.get_column("contains_tryptophan").to_list() == [False, None, False]
