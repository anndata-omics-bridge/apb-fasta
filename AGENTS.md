# APB FASTA — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --group dev --group docs` |
| Format | `.venv/bin/ruff format src tests && .venv/bin/ruff check --fix src tests` |
| Lint | `.venv/bin/ruff check src tests` |
| Architecture | `.venv/bin/lint-imports` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Docs | `uv run --group docs zensical build --clean --strict` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Full gate | `make check` |

## Architecture

- `api.py` owns the public protein-bound lifecycle; each operation accepts its APB2 result explicitly.
- `integration.py` is the only module that translates or updates APB2 result values.
- `calculation/` accepts and returns ordinary Polars values and imports no APB2 or storage framework.
- `cli.py` only composes APB2 result I/O and the public annotation lifecycle from `apb_fasta.api`.
- `protein_fasta` owns FASTA reading and header interpretation; Prozor owns peptide matching; APB2 owns result persistence.

## Code conventions

- Fully annotate every function and method in `src/` and `tests`.
- Use strict Pyright and Ruff with the configured 100-character line length.
- Keep `__init__.py` empty and import from defining modules inside this package. Other anndata_bridge packages import this one only from `apb_fasta.api`, and it imports them only from theirs. The CLI imports this package only from `apb_fasta.api` too.
- Use Google-style docstrings for public APIs.
- Keep `docs/` current with behavioral changes; `make check` builds it with `--strict`.
- Preserve unrelated worktree changes and add focused tests with behavioral changes.
