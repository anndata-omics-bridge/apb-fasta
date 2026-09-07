# Verify peptides and annotate proteins

Both operations need the same two inputs: an APB2 result and a parsed protein database.

## Load the inputs once

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels, write_parsed_levels
from apb_fasta.annotation import FastaAnnotationParser
from protein_fasta.frame import ProteinDatabase, refseq, uniprotkb

proteins = ProteinDatabase(uniprotkb, refseq).parse(
    (Path("human.fasta"), Path("contaminants.fasta"))
)
anndata = read_parsed_levels(Path("input.h5mu"))
parser = FastaAnnotationParser(anndata, proteins)
```

The protein frame is parsed once and reused by every operation. Header formats are supplied in priority order.

## Verify peptides only

```python
verified = parser.verify_peptides()
print(verified.reports.peptide_levels)
write_parsed_levels(verified.parsed, Path("verified.h5mu"))
```

Verification covers every peptide-derived level that carries APB2's canonical `ProForma_peptide` column. A result with no such level is an error, not a silent no-op.

## Merge protein annotations only

```python
annotated = parser.merge_annotations()
print(annotated.reports.protein_groups)
write_parsed_levels(annotated.parsed, Path("annotated.h5mu"))
```

Annotation reads the protein axis and the FASTA-accession column role APB2 persisted for it. A result with no protein level is an error.

## Apply both

```python
complete = parser.run()
write_parsed_levels(complete.parsed, Path("complete.h5mu"))
```

`run` chains the two operations in memory and returns both reports. Chaining them by hand is equivalent:

```python
verified = parser.verify_peptides()
complete = FastaAnnotationParser(verified.parsed, proteins).merge_annotations()
```

## Configuration

```python
from apb_fasta.configuration import FastaAnnotationParameters

parameters = FastaAnnotationParameters(
    protein_group_separator=";",
    matcher_backend="auto",
    il_equivalent=False,
)
parser = FastaAnnotationParser(anndata, proteins, parameters=parameters)
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `protein_group_separator` | `";"` | Separator splitting a protein group into members; must not be empty |
| `matcher_backend` | `"auto"` | Prozor Aho-Corasick backend: `auto`, `ahocorapy`, or `ahocorasick_rs` |
| `il_equivalent` | `False` | Treat leucine and isoleucine as equivalent when matching |

## Failures are explicit

Every method returns an immutable `FastaAnnotationResult`. Missing levels, missing protein-frame columns, and existing output names that would be overwritten all raise `FastaAnnotationError` before any replacement result is constructed. Nothing is written on a partial application.
