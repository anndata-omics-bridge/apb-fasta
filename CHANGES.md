# Changes

- 2026-10-02: Added `add_peptide_properties(parsed)`, which writes a feature-aligned `varm["peptide_properties"]` with protein_fasta's sequence-derived peptide properties (length, mass, pI, hydrophobicity, instability and Boman indices, charge, predicted reversed-phase retention time, missed cleavages, motif flags) on every peptide-derived level, and records the protein_fasta version under `fasta.provenance.peptide_properties`. It is not part of `annotate()` or the CLI.

- 2026-10-02: `varm["fasta_validation"]` gains `fasta_matching_organisms` (sorted, `;`-joined organism mnemonics of the matching proteins) and `fasta_matches_contaminant` (a matching protein has `is_contaminant`), taken from protein_fasta's parsed columns.

- 2026-09-19: **Breaking:** Replaced dataset-bound `FastaAnnotationParser(parsed, proteins)` with reusable `FastaAnnotator(proteins)`. `verify_peptides(parsed)`, `merge_annotations(parsed)`, and `annotate(parsed)` now accept canonical APB2 values explicitly; the Python `run()` method and `apb_fasta.annotation` module were removed. CLI command names and scientific results are unchanged.

- 2026-09-18: **Breaking:** FASTA metadata schema 2 moves root database sources and operation settings into `fasta.provenance`, grouped by operation. Per-level validation summaries and aligned tables are unchanged. APB2 stores the contributions on their owning objects and combines them only for standalone H5AD.

- 2026-09-17: FASTA configuration and source provenance are stored once in shared metadata; per-level peptide sections retain only their level-specific coverage report.

- 2026-09-03: Introduced the APB2-bound FASTA lifecycle. The explicit `verify_peptides()` and `merge_annotations()` operations can run independently, compose in memory through `run()`, and are exposed as matching CLI commands. The package accepts the configured Polars frame produced by `protein_fasta`, validates canonical stripped peptide sequences with Prozor, expands all reported protein-group members into a lossless long annotation table, and returns a copied `ParsedLevels` result with operation-specific reports and provenance. Added H5AD, H5MU, Parquet, and DuckDB coverage. Prozor protein inference remains a later opt-in operation.
