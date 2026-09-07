.PHONY: format format-check lint typecheck test check build

format:
	.venv/bin/ruff format src tests
	.venv/bin/ruff check --fix src tests

format-check:
	.venv/bin/ruff format --check src tests

lint:
	.venv/bin/ruff check src tests
	.venv/bin/lint-imports

typecheck:
	.venv/bin/pyright

test:
	.venv/bin/pytest --cov --cov-branch

build:
	uv build
	.venv/bin/twine check dist/*

check:
	uv lock --check
	$(MAKE) format-check lint typecheck test build
	.venv/bin/deptry .
