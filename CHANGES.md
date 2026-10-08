# Changes

## 2026-10-08

- Peptide coverage separates decoys: `PeptideCoverage` gains `decoy_feature_count`, the rows apb2 marks `apb_Decoy`, and `matched_feature_count` and `unmatched_feature_count` now count targets only, so a decoy no longer reports as an unmatched peptide. Peptide verification requires `apb_Decoy` on every level it checks; the CLI reports `decoys=`.

## 2026-10-05

- `uv.lock` is no longer committed: `make sync` and CI resolve the environment from `pyproject.toml`, `make check` drops `uv lock --check`, CI caches by `pyproject.toml`, and the dev group pins `ruff==0.16.10` and `pyright==1.1.414` so lint and type results stay stable.
- README shows the Zenodo DOI badge, linking the concept DOI that resolves to the latest archived release, and the PyPI version badge.
- 0.1.1: `CITATION.cff` records the author's ORCID, so the Zenodo archive of each GitHub release carries complete citation metadata and a DOI. CI checks out APB2 0.1.1, protein-fasta 0.3.1 and Prozor 0.1.1.
- PyPI release setup: `.github/workflows/publish.yml` builds and checks the distributions, then publishes to PyPI through trusted publishing for a published GitHub release tagged `v<version>`; a manual run only builds and checks. README links point at GitHub or the documentation site so they resolve on PyPI, README gains an installation section, and `pyproject.toml` adds keywords and classifiers.

## 2026-10-04

- `FastaAnnotator`, its methods and `add_peptide_properties` drop the `/` and `*` signature markers; every call that worked before still works.
- `FastaAnnotator.read(fasta, formats, parameters)` reads FASTA files, or the protein-fasta database Parquet, so callers no longer build the protein table with protein_fasta; `FastaAnnotator.proteins` returns the bound table. The CLI uses `read`, and an unknown `--formats` name now fails with protein_fasta's own message.

## 2026-10-03

- Peptide matching searches 10,000 protein sequences per Aho-Corasick call through prozor and builds no per-protein record objects; about 2.7 times faster on ProteoBench's 2.84 M-entry entrapment FASTA. FASTA_PATHS may instead be one protein-fasta database Parquet file (`protein-fasta database`), read about 70 times faster than parsing.
- Read FASTA-accession evidence from VarFinal.roles.

- 2026-10-02: Added `add_peptide_properties(parsed)`, which writes a feature-aligned `varm["peptide_properties"]` with protein_fasta's sequence-derived peptide properties (length, mass, pI, hydrophobicity, instability and Boman indices, charge, predicted reversed-phase retention time, missed cleavages, motif flags) on every peptide-derived level, and records the protein_fasta version under `fasta.provenance.peptide_properties`. It is not part of `annotate()` or the CLI.

- 2026-10-02: `varm["fasta_validation"]` gains `fasta_matching_organisms` (sorted, `;`-joined organism mnemonics of the matching proteins) and `fasta_matches_contaminant` (a matching protein has `is_contaminant`), taken from protein_fasta's parsed columns.

- 2026-09-19: **Breaking:** Replaced dataset-bound `FastaAnnotationParser(parsed, proteins)` with reusable `FastaAnnotator(proteins)`. `verify_peptides(parsed)`, `merge_annotations(parsed)`, and `annotate(parsed)` now accept canonical APB2 values explicitly; the Python `run()` method and `apb_fasta.annotation` module were removed. CLI command names and scientific results are unchanged.

- 2026-09-18: **Breaking:** FASTA metadata schema 2 moves root database sources and operation settings into `fasta.provenance`, grouped by operation. Per-level validation summaries and aligned tables are unchanged. APB2 stores the contributions on their owning objects and combines them only for standalone H5AD.

- 2026-09-17: FASTA configuration and source provenance are stored once in shared metadata; per-level peptide sections retain only their level-specific coverage report.

- 2026-09-03: Introduced the APB2-bound FASTA lifecycle. The explicit `verify_peptides()` and `merge_annotations()` operations can run independently, compose in memory through `run()`, and are exposed as matching CLI commands. The package accepts the configured Polars frame produced by `protein_fasta`, validates canonical stripped peptide sequences with Prozor, expands all reported protein-group members into a lossless long annotation table, and returns a copied `ParsedLevels` result with operation-specific reports and provenance. Added H5AD, H5MU, Parquet, and DuckDB coverage. Prozor protein inference remains a later opt-in operation.
