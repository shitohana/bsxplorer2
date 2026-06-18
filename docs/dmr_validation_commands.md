# DMR Validation Commands

DMR validation business logic lives inside package layers:

```text
python/src/dmr_validation_framework/checks/
python/src/dmr_validation_framework/workflows/
python/src/dmr_validation_framework/importers/
python/src/dmr_validation_framework/models/
python/src/dmr_validation_framework/reports/
```

The repository-level `scripts/` directory is intentionally empty.

Use the DMR validation framework entry point for validation checks:

```bash
PYTHONPATH=python/src python3 -m dmr_validation_framework.cli list
PYTHONPATH=python/src python3 -m dmr_validation_framework.cli profiles
PYTHONPATH=python/src python3 -m dmr_validation_framework.cli run-all --profile core --out-dir outputs/validation_audit
```

Profiles:

- `core`: default DMR validation: input audit, coverage balance, effect robustness, caller agreement, GLM/GLMM agreement, model status.
- `extended`: `core` plus optional checks.
- `extended_only`: optional checks only.
- `case_study`: reserved for clean case-study dispatchers; generic run-all does not include Oryza-specific workflows.

New DMR validation logic should go under:

```text
python/src/dmr_validation_framework/
```

Dataset-specific workflows should go under:

```text
case_studies/
```

The Oryza GSE202715 downstream workflow has been moved to:

```text
case_studies/oryza_gse202715/functional_downstream/
```
