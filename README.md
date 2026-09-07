# apb-fasta

FASTA verification and protein annotation for APB2 results.

**[Online documentation](https://anndata-omics-bridge.github.io/apb-fasta/)** or its [source index](docs/index.md).

The first release exposes two independent operations:

- verify that APB2's modification-stripped peptide sequences occur in the supplied FASTA database;
- merge FASTA annotations for every reported protein-group member without collapsing the group to its leading accession.

Protein inference with Prozor is the planned third operation. It will remain explicit and opt-in.

## Python API

Load the APB2 result and parse the protein database once:

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

Verify peptides only:

```python
verified = parser.verify_peptides()
print(verified.reports.peptide_levels)
write_parsed_levels(verified.parsed, Path("verified.h5mu"))
```

Merge protein annotations only:

```python
annotated = parser.merge_annotations()
print(annotated.reports.protein_groups)
write_parsed_levels(annotated.parsed, Path("annotated.h5mu"))
```

Apply both operations in memory:

```python
complete = parser.run()
write_parsed_levels(complete.parsed, Path("complete.h5mu"))
```

The operations can also be chained explicitly without an intermediate file:

```python
verified = parser.verify_peptides()
complete = FastaAnnotationParser(verified.parsed, proteins).merge_annotations()
```

Every method returns an immutable `FastaAnnotationResult` containing the replacement APB2 `ParsedLevels` value and typed operation reports. `apb_fasta` neither opens result files nor receives raw AnnData or MuData objects.

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

Peptide verification adds feature-aligned `varm["fasta_validation"]` tables. Protein annotation adds protein-aligned `varm["fasta"]`, the lossless `fasta_protein_group_members` annotation table, and its directed relation to the protein axis.
