# APB FASTA — agent rules

The closest `AGENTS.md` wins. Explicit user instructions override this file.

## Verified commands

| Task | Command |
| --- | --- |
| Synchronize | `uv sync --frozen --group dev --group docs` |
| Format | `.venv/bin/ruff format src tests && .venv/bin/ruff check --fix src tests` |
| Lint | `.venv/bin/ruff check src tests` |
| Architecture | `.venv/bin/lint-imports` |
| Typecheck | `.venv/bin/pyright` |
| Dependencies | `.venv/bin/deptry .` |
| Tests | `.venv/bin/pytest --cov --cov-branch` |
| Docs | `uv run --frozen --group docs zensical build --clean --strict` |
| Build | `uv build && .venv/bin/twine check dist/*` |
| Full gate | `make check` |

## Architecture

- `annotation.py` owns the public dataset-and-protein-bound lifecycle.
- `integration.py` is the only module that translates or updates APB2 result values.
- `calculation/` accepts and returns ordinary Polars values and imports no APB2 or storage framework.
- `cli.py` only composes `protein_fasta`, APB2 result I/O, and the public annotation lifecycle.
- `protein_fasta` owns FASTA reading and header interpretation; Prozor owns peptide matching; APB2 owns result persistence.

## Code conventions

- Fully annotate every function and method in `src/` and `tests`.
- Use strict Pyright and Ruff with the configured 100-character line length.
- Keep `__init__.py` empty and import public objects from defining modules.
- Use Google-style docstrings for public APIs.
- Keep `docs/` current with behavioral changes; `make check` builds it with `--strict`.
- Preserve unrelated worktree changes and add focused tests with behavioral changes.
