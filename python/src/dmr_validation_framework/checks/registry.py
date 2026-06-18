"""Registry for lightweight DMR validation checks."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CheckSpec:
    name: str
    module: str
    layer: str
    outputs: tuple[str, ...]
    description: str
    profiles: tuple[str, ...] = ("extended_only",)

    @property
    def entrypoint(self) -> str:
        return f"{self.module}:main"


CHECKS: tuple[CheckSpec, ...] = (
    CheckSpec("data_inventory", "dmr_validation_framework.checks.input_audit", "input_audit", ("input_inventory.tsv",), "Inventory available DMR validation inputs and required columns.", ("core",)),
    CheckSpec("coverage_set_balance", "dmr_validation_framework.checks.coverage_set_balance", "coverage_qc", ("coverage_set_summary.tsv",), "Audit per-sample vs common-CpG coverage-set balance.", ("core",)),
    CheckSpec("delta_weighting_sensitivity", "dmr_validation_framework.checks.delta_weighting", "effect_robustness", ("delta_weighting_summary.tsv",), "Compare replicate-weighted and pooled-read delta estimates with LOO/bootstrap diagnostics.", ("core",)),
    CheckSpec("delta_bootstrap_ci", "dmr_validation_framework.checks.delta_bootstrap", "effect_robustness", ("delta_bootstrap_summary.tsv",), "Bootstrap confidence intervals for replicate-level delta methylation.", ("core",)),
    CheckSpec("direction_agreement", "dmr_validation_framework.checks.direction_agreement", "caller_agreement", ("direction_agreement.tsv",), "Evaluate matched-caller delta direction agreement.", ("core",)),
    CheckSpec("overlap_threshold_sensitivity", "dmr_validation_framework.checks.overlap_thresholds", "caller_agreement", ("overlap_threshold_sensitivity.tsv",), "Evaluate caller overlap sensitivity across reciprocal-overlap thresholds.", ("core",)),
    CheckSpec("glm_glmm_validation_status", "dmr_validation_framework.checks.glm_glmm_status", "model_validation", ("glm_glmm_status_summary.tsv",), "Summarize GLM vs GLMM significance and status agreement.", ("core",)),
    CheckSpec("model_status_audit", "dmr_validation_framework.checks.model_status", "model_validation", ("model_status_audit.tsv",), "Audit model convergence, fallback, and failure status flags.", ("core",)),
    CheckSpec("projection_sensitivity", "dmr_validation_framework.checks.projection_sensitivity", "metagene_qc", ("projection_sensitivity_summary.tsv",), "Compare DMR center and interval projection strategies.", ("extended_only",)),
    CheckSpec("density_denominator_sensitivity", "dmr_validation_framework.checks.density_denominator", "metagene_qc", ("density_denominator_sensitivity.tsv",), "Compare DMR occupancy density denominator choices.", ("extended_only",)),
    CheckSpec("random_control_occupancy", "dmr_validation_framework.checks.random_control_occupancy", "negative_control", ("random_control_occupancy_summary.tsv",), "Compare observed DMR occupancy with random-control regions.", ("extended_only",)),
)

PROFILE_DESCRIPTIONS = {
    "core": "Default DMR validation: input audit, coverage balance, effect robustness, caller agreement, GLM/GLMM agreement, model status.",
    "extended": "Core profile plus all optional validation checks.",
    "extended_only": "Only optional checks not included in the default core profile.",
    "case_study": "Reserved for clean case-study dispatchers; generic run-all does not include dataset-specific workflows.",
}


def check_names() -> list[str]:
    return [check.name for check in CHECKS]


def profile_names() -> list[str]:
    return list(PROFILE_DESCRIPTIONS)


def checks_for_profile(profile: str) -> tuple[CheckSpec, ...]:
    if profile == "extended":
        wanted = {"core", "extended_only"}
    elif profile in PROFILE_DESCRIPTIONS:
        wanted = {profile}
    else:
        raise KeyError(f"Unknown DMR validation profile: {profile}")
    return tuple(check for check in CHECKS if set(check.profiles) & wanted)


def get_check(name: str) -> CheckSpec:
    for check in CHECKS:
        if check.name == name:
            return check
    raise KeyError(f"Unknown DMR validation check: {name}")
