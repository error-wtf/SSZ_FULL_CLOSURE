.PHONY: test audit full direct clean

test:
	pytest -q

audit:
	python ssz_p5_full_closure_auditor.py --data-dir . --quick

full:
	python ssz_p5_full_closure_auditor.py --data-dir . --full

direct:
	python ssz_p5_full_closure_auditor.py --data-dir . --full --require-direct-krgm

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache .mypy_cache build dist *.egg-info src/*.egg-info

full-pipeline:
	python ssz_p5_full_pipeline.py

full-pipeline-strict:
	python ssz_p5_full_pipeline.py --strict
