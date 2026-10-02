"""Feature-aligned sequence properties of one peptide-derived level."""

from __future__ import annotations

import polars as pl
from protein_fasta.peptide_frame import peptide_property_frame


def peptide_level_properties(frame: pl.DataFrame, sequence_column: str, /) -> pl.DataFrame:
    """Return one property row per feature of ``frame``, in feature order.

    Each distinct sequence is computed once; a feature without a sequence has null properties.
    """
    features = frame.select(
        pl.col(sequence_column)
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase()
        .alias("sequence")
    ).select(pl.when(pl.col("sequence") != "").then(pl.col("sequence")).alias("sequence"))
    sequences = features.get_column("sequence").drop_nulls().unique(maintain_order=True)
    properties = peptide_property_frame(sequences.to_list())
    return features.join(properties, on="sequence", how="left", maintain_order="left").drop(
        "sequence"
    )
