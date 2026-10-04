"""Feature-aligned peptide-property calculation."""

from __future__ import annotations

from dataclasses import fields

import polars as pl
import pytest
from protein_fasta.api import PeptideProperties

from apb_fasta.calculation.peptide_properties import peptide_level_properties


def test_properties_align_with_features_including_duplicates_and_missing() -> None:
    frame = pl.DataFrame(
        {"ProForma_peptide": [" kqpwwr ", None, "", "EAAAMGPTK", "KQPWWR", "RSGASILQAGCUG"]}
    )

    properties = peptide_level_properties(frame, "ProForma_peptide")

    assert properties.height == frame.height
    assert properties.columns == [field.name for field in fields(PeptideProperties)]
    assert properties.get_column("length").to_list() == [6, None, None, 9, 6, None]
    assert properties.get_column("c_terminal_residue").to_list() == [
        "R",
        None,
        None,
        "K",
        "R",
        None,
    ]


def test_modified_sequence_is_rejected() -> None:
    frame = pl.DataFrame({"ProForma_peptide": ["PEPT[Phospho]IDE"]})

    with pytest.raises(ValueError, match="upper-case letters"):
        peptide_level_properties(frame, "ProForma_peptide")
