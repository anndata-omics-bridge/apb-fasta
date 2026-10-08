# Verify peptides and annotate proteins

Both operations need the same two inputs: an APB2 result and the FASTA files.

## Load the inputs once

```python
from pathlib import Path

from apb2.api import read_parsed_levels, write_parsed_levels
from apb_fasta.api import FastaAnnotator

annotator = FastaAnnotator.read((Path("human.fasta"), Path("contaminants.fasta")))
parsed = read_parsed_levels(Path("input.h5mu"))
```

`FastaAnnotator.read()` parses the FASTA files once, with protein_fasta's `uniprotkb` and `refseq` header formats unless `formats` names others in priority order, and reuses the table for every operation. It also accepts the Parquet database `protein-fasta database` writes; `proteins` returns the bound table. `FastaAnnotator(proteins)` binds a table that is already parsed.

## Verify peptides only

```python
verified = annotator.verify_peptides(parsed)
print(verified.reports.peptide_levels)
write_parsed_levels(verified.parsed, Path("verified.h5mu"))
```

Verification covers every peptide-derived level that carries APB2's canonical `ProForma_peptide` column. A result with no such level is an error, not a silent no-op.

## Merge protein annotations only

```python
annotated = annotator.merge_annotations(parsed)
print(annotated.reports.protein_groups)
write_parsed_levels(annotated.parsed, Path("annotated.h5mu"))
```

Annotation reads the protein axis and the FASTA-accession column role APB2 persisted for it. A result with no protein level is an error.

## Apply both

```python
complete = annotator.annotate(parsed)
write_parsed_levels(complete.parsed, Path("complete.h5mu"))
```

`annotate` chains the two operations in memory and returns both reports. Chaining them by hand is equivalent:

```python
verified = annotator.verify_peptides(parsed)
complete = annotator.merge_annotations(verified.parsed)
```

## Configuration

```python
from apb_fasta.api import FastaAnnotationParameters

parameters = FastaAnnotationParameters(
    protein_group_separator=";",
    matcher_backend="auto",
)
annotator = FastaAnnotator(proteins, parameters=parameters)
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `protein_group_separator` | `";"` | Separator splitting a protein group into members; must not be empty |
| `matcher_backend` | `"auto"` | Prozor Aho-Corasick backend: `auto`, `ahocorapy`, or `ahocorasick_rs` |

## Failures are explicit

Every method returns an immutable `FastaAnnotationResult`. Missing levels, missing protein-frame columns, and existing output names that would be overwritten all raise `FastaAnnotationError` before any replacement result is constructed. Nothing is written on a partial application.
