"""Public FASTA annotation lifecycle and APB2 persistence."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from apb2.result_facade import (
    FinalLayerTable,
    JsonValue,
    ObsFinal,
    ParsedLevel,
    ParsedLevels,
    VarFinal,
    read_parsed_levels,
    write_parsed_levels,
)

from apb_fasta.annotation import FastaAnnotationParser
from apb_fasta.cli import app
from apb_fasta.errors import FastaAnnotationError


def _level(
    name: str,
    var: pl.DataFrame,
    key_columns: tuple[str, ...],
    column_roles: dict[str, str],
) -> ParsedLevel:
    layer_name = "Intensity"
    roles: dict[str, JsonValue] = {}
    for key, value in column_roles.items():
        roles[key] = value
    uns: dict[str, JsonValue] = {
        "quantification_level": name,
        "column_roles": roles,
        "matrix_values_projected": True,
    }
    return ParsedLevel(
        obs=ObsFinal(frame=pl.DataFrame({"Run": ["run1"]}), key_columns=("Run",)),
        var=VarFinal(frame=var, key_columns=key_columns),
        primary_layer_name=layer_name,
        uns=uns,
        layers={
            layer_name: FinalLayerTable(
                layer_name=layer_name,
                var_key_columns=key_columns,
                values=pl.DataFrame(
                    {
                        **{key: var.get_column(key) for key in key_columns},
                        "obs_0": [float(index + 1) for index in range(var.height)],
                    }
                ),
            )
        },
        obsm={},
        varm={},
        obsp={},
        varp={},
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


def test_run_verifies_peptides_and_merges_annotations_without_mutating_input() -> None:
    source = _parsed()

    result = FastaAnnotationParser(source, _proteins()).run()

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
    assert set(metadata) == {
        "schema_version",
        "peptide_verification",
        "protein_annotation",
    }


def test_operations_run_independently() -> None:
    parser = FastaAnnotationParser(_parsed(), _proteins())

    verified = parser.verify_peptides()
    annotated = parser.merge_annotations()

    assert "fasta_validation" in verified.parsed.levels["peptide"].varm
    assert verified.parsed.levels["protein"].varm == {}
    assert verified.parsed.annotation_tables == {}
    assert verified.reports.protein_groups is None
    assert annotated.parsed.levels["peptide"].varm == {}
    assert "fasta" in annotated.parsed.levels["protein"].varm
    assert annotated.reports.peptide_levels == {}
    assert annotated.reports.protein_groups is not None


def test_operations_compose_in_memory() -> None:
    verified = FastaAnnotationParser(_parsed(), _proteins()).verify_peptides()

    complete = FastaAnnotationParser(verified.parsed, _proteins()).merge_annotations()

    assert "fasta_validation" in complete.parsed.levels["peptide"].varm
    assert "fasta" in complete.parsed.levels["protein"].varm
    metadata = complete.parsed.metadata["fasta"]
    assert isinstance(metadata, dict)
    assert set(metadata) == {
        "schema_version",
        "peptide_verification",
        "protein_annotation",
    }


def test_operations_compose_in_reverse_order() -> None:
    annotated = FastaAnnotationParser(_parsed(), _proteins()).merge_annotations()

    complete = FastaAnnotationParser(annotated.parsed, _proteins()).verify_peptides()

    assert "fasta_validation" in complete.parsed.levels["peptide"].varm
    assert "fasta" in complete.parsed.levels["protein"].varm


@pytest.mark.parametrize("suffix", [".h5mu", ".parquet", ".duckdb"])
def test_annotated_result_round_trips_through_multilevel_formats(
    suffix: str,
    tmp_path: Path,
) -> None:
    result = FastaAnnotationParser(_parsed(), _proteins()).run().parsed
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
    result = FastaAnnotationParser(peptide_only, _proteins()).verify_peptides().parsed
    target = tmp_path / "annotated.h5ad"

    write_parsed_levels(result, target)
    restored = read_parsed_levels(target)

    assert restored.levels["peptide"].varm["fasta_validation"].height == 3
    assert restored.annotation_tables == {}
    assert restored.feature_relations == {}


def test_parser_rejects_owned_output_collisions() -> None:
    parsed = _parsed()
    parsed.levels["peptide"].varm["fasta_validation"] = pl.DataFrame({"old": [1, 2, 3]})

    with pytest.raises(FastaAnnotationError, match="already contains"):
        FastaAnnotationParser(parsed, _proteins()).verify_peptides()


def test_parser_rejects_repeating_an_applied_operation() -> None:
    verified = FastaAnnotationParser(_parsed(), _proteins()).verify_peptides()

    with pytest.raises(FastaAnnotationError, match="peptide_verification"):
        FastaAnnotationParser(verified.parsed, _proteins()).verify_peptides()


def test_parser_requires_protein_fasta_database_schema() -> None:
    with pytest.raises(FastaAnnotationError, match="missing required columns"):
        FastaAnnotationParser(
            _parsed(),
            _proteins().drop("fasta_source_checksum"),
        ).verify_peptides()


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
