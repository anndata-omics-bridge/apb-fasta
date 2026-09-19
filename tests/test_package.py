"""Package surface tests."""

from __future__ import annotations

import ast
from pathlib import Path


def test_package_imports() -> None:
    import apb_fasta.api
    import apb_fasta.cli
    import apb_fasta.configuration

    assert apb_fasta.api.FastaAnnotator
    assert apb_fasta.cli.app
    assert apb_fasta.configuration.FastaAnnotationParameters


def test_public_api_is_an_in_memory_annotation_boundary() -> None:
    path = Path(__file__).parents[1] / "src/apb_fasta/api.py"
    document = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imported = {
        (node.module, name.name)
        for node in ast.walk(document)
        if isinstance(node, ast.ImportFrom) and node.module is not None
        for name in node.names
    }

    assert {item for item in imported if item[0].startswith("apb2")} == {
        ("apb2.api", "ParsedLevels")
    }
    assert not any(
        module == "pathlib" or module.startswith("protein_fasta") for module, _name in imported
    )
    assert not {"read_parsed_levels", "write_parsed_levels"}.intersection(
        name for _module, name in imported
    )
