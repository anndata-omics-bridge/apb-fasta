# apb-fasta

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23151458.svg)](https://doi.org/10.5281/zenodo.23151458)
[![PyPI](https://img.shields.io/pypi/v/apb-fasta.svg)](https://pypi.org/project/apb-fasta/)

FASTA verification and protein annotation for APB2 results.

**[Online documentation](https://anndata-omics-bridge.github.io/apb-fasta/)** or its [source index](https://github.com/anndata-omics-bridge/apb-fasta/blob/main/docs/index.md).

The first release exposes two independent operations:

- verify that APB2's modification-stripped peptide sequences occur in the supplied FASTA database;
- merge FASTA annotations for every reported protein-group member without collapsing the group to its leading accession.

Protein inference with Prozor is the planned third operation. It will remain explicit and opt-in.

## Installation

APB FASTA requires Python 3.13 or later.

```bash
pip install apb-fasta
```

## Python API

Read the FASTA files once and load the APB2 result:

```python
from pathlib import Path

from apb2.api import read_parsed_levels, write_parsed_levels
from apb_fasta.api import FastaAnnotator

annotator = FastaAnnotator.read((Path("human.fasta"), Path("contaminants.fasta")))
parsed = read_parsed_levels(Path("input.h5mu"))
```

Verify peptides only:

```python
verified = annotator.verify_peptides(parsed)
print(verified.reports.peptide_levels)
write_parsed_levels(verified.parsed, Path("verified.h5mu"))
```

Merge protein annotations only:

```python
annotated = annotator.merge_annotations(parsed)
print(annotated.reports.protein_groups)
write_parsed_levels(annotated.parsed, Path("annotated.h5mu"))
```

Apply both operations in memory:

```python
complete = annotator.annotate(parsed)
write_parsed_levels(complete.parsed, Path("complete.h5mu"))
```

The operations can also be chained explicitly without an intermediate file:

```python
verified = annotator.verify_peptides(parsed)
complete = annotator.merge_annotations(verified.parsed)
```

The annotator validates and binds the reusable protein frame once. Every method accepts one canonical `ParsedLevels` value and returns an immutable `FastaAnnotationResult` containing its replacement plus typed operation reports. `apb_fasta` neither opens result files nor receives raw AnnData or MuData objects.

## CLI

Verify peptides only:

```bash
apb-fasta verify-peptides input.h5mu human.fasta contaminants.fasta --output verified.h5mu
```

Merge protein annotations only:

```bash
apb-fasta merge-annotations input.h5mu human.fasta contaminants.fasta --output annotated.h5mu
```

Apply both operations together:

```bash
apb-fasta run input.h5mu human.fasta contaminants.fasta --output complete.h5mu
```

Peptide verification adds feature-aligned `varm["fasta_validation"]` tables. Its coverage counts matched and unmatched targets apart from decoys, the rows apb2 marks `apb_Decoy`, which every checked level must carry. Protein annotation adds protein-aligned `varm["fasta"]`, the lossless `fasta_protein_group_members` annotation table, and its directed relation to the protein axis.
