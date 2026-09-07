# Command-line interface

```bash
apb-fasta --help
```

Three commands, one per operation. All take an APB2 result, one or more FASTA files, and a mandatory `--output`.

## `verify-peptides`

```bash
apb-fasta verify-peptides input.h5mu human.fasta contaminants.fasta --output verified.h5mu
```

## `merge-annotations`

```bash
apb-fasta merge-annotations input.h5mu human.fasta contaminants.fasta --output annotated.h5mu
```

## `run`

```bash
apb-fasta run input.h5mu human.fasta contaminants.fasta --output complete.h5mu
```

Applies both operations in one in-memory pass and reports both coverage summaries.

## Options

| Option | Default | Commands |
| --- | --- | --- |
| `--output` | required | all |
| `--formats` | `uniprotkb refseq` | all |
| `--protein-group-separator` | `;` | all |
| `--backend` | `auto` | `verify-peptides`, `run` |
| `--il-equivalent` | off | `verify-peptides`, `run` |

`--formats` names the `protein_fasta` header formats to try, in priority order. `--backend` selects the Prozor Aho-Corasick implementation.

## Exit codes

`0` on success. `1` on an expected failure — a missing file, a missing level, a protein frame without the required columns, or an output name that already exists. The reason is logged; nothing is written.
