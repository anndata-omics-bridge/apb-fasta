# APB FASTA

FASTA verification and protein annotation for [APB2](https://anndata-omics-bridge.github.io/apb2/) results.

`apb-fasta` takes one APB2 result and one parsed protein database and returns a replacement result. It never opens result files itself, and it never receives raw AnnData or MuData objects — [`protein_fasta`](https://anndata-omics-bridge.github.io/protein-fasta/) owns FASTA reading and header interpretation, [Prozor](https://anndata-omics-bridge.github.io/prozor/) owns peptide matching, and APB2 owns result persistence.

## Two independent operations

| Operation | What it does |
| --- | --- |
| `verify_peptides` | Confirms every modification-stripped peptide sequence occurs in the supplied FASTA database |
| `merge_annotations` | Merges FASTA annotations for every reported protein-group member, without collapsing the group to its leading accession |

Each operation is independently callable and independently persistable. `run` applies both in memory without an intermediate file.

Protein inference with Prozor is the planned third operation. It will remain explicit and opt-in.

## Install

`apb-fasta` is developed alongside its siblings in one workspace:

```bash
uv sync --frozen --group dev --group docs
```

## Next

- [Verify peptides and annotate proteins](workflow.md) — the end-to-end guide
- [Result layout](results.md) — what the operations write into the result
- [Command-line interface](cli.md) and [Python API](api.md)
- [Architecture](architecture.md) — module boundaries and the import direction
