from dmr_validation_framework.checks import checks_for_profile


def names(profile: str) -> set[str]:
    return {check.name for check in checks_for_profile(profile)}


def test_core_profile_is_minimal_default_validation_set():
    assert names("core") == {
        "data_inventory",
        "coverage_set_balance",
        "delta_weighting_sensitivity",
        "delta_bootstrap_ci",
        "direction_agreement",
        "overlap_threshold_sensitivity",
        "glm_glmm_validation_status",
        "model_status_audit",
    }


def test_extended_profile_is_core_plus_optional():
    assert names("extended") == names("core") | names("extended_only")


def test_extended_only_profile_contains_optional_checks():
    assert names("extended_only") == {
        "density_denominator_sensitivity",
        "projection_sensitivity",
        "random_control_occupancy",
    }
    assert "density_denominator_sensitivity" not in names("core")


def test_case_study_profile_does_not_run_oryza_specific_workflows():
    assert names("case_study") == set()
