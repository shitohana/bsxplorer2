#!/usr/bin/env python3
"""Run all lightweight critical validation checks for DMR/metagene outputs."""

from __future__ import annotations

import argparse
import importlib
import io
from pathlib import Path
from contextlib import redirect_stderr, redirect_stdout

from dmr_validation_framework.checks import CHECKS as REGISTERED_CHECKS  # noqa: E402
from dmr_validation_framework.checks import PROFILE_DESCRIPTIONS, checks_for_profile, profile_names  # noqa: E402
from dmr_validation_framework.core.io import ensure_out_dir, read_table, write_tsv  # noqa: E402

EXPECTED_OUTPUTS = {check.name: list(check.outputs) for check in REGISTERED_CHECKS}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--profile", choices=profile_names(), default="core")
    parser.add_argument("--list-profiles", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def infer_status(out_dir: Path, check_name: str, returncode: int) -> tuple[str, str]:
    if returncode != 0:
        return "FAIL", f"script exited with code {returncode}"
    expected = EXPECTED_OUTPUTS[check_name]
    paths = [out_dir / p for p in expected]
    if not all(p.exists() for p in paths):
        return "FAIL", "expected output file missing"
    text = "\n".join(p.read_text(encoding="utf-8", errors="replace")[:4000] for p in paths)
    if "SKIPPED" in text:
        return "SKIPPED", "check skipped with explicit reason"
    if "WARN" in text or "missing_status" in text:
        return "WARN", "check completed with warnings"
    return "PASS", "check completed"


def write_acceptance(out_dir: Path) -> None:
    text = """# Acceptance criteria for critical DMR/metagene validation audit

- Все проверки имеют статус PASS/WARN или SKIPPED с явной причиной.
- Нет молчаливых пропусков: каждый ожидаемый output создан либо содержит SKIPPED.
- Все TSV outputs имеют заголовки.
- Все runtime outputs записываются только в `outputs/validation_audit/`.
- Скрипты не запускают FASTQ/Bismark/raw alignment/DSS/methylKit/dmrseq/metilene pipelines.
- Отчёт не содержит claims о genome-wide DMR caller.
- GLMM описан только как confirmatory validation для выбранных регионов.
- q-values DSS, methylKit и внутренних моделей не объявлены статистически эквивалентными.
- coverage_set_balance выдаёт PASS/WARN/SKIPPED явно.
- delta_weighting_sensitivity выдаёт PASS/WARN/SKIPPED явно.
- delta_bootstrap_ci выдаёт PASS/WARN/SKIPPED явно.
- bootstrap CI записан как диагностическая неопределённость effect size, а не p-value или genome-wide DMR significance.
- При малом числе биологических реплик присутствует small-n warning, если применимо.
- random_control output содержит `empirical_q_BH` или явный SKIPPED/WARN.
- random_control output содержит `global_profile_p` или явный SKIPPED/WARN.
- direction agreement содержит `delta_threshold`.
- overlap sensitivity содержит `matching_policy`.
- overlap filtering uses strict reciprocal overlap `min(O/L_a, O/L_b)`.
- Legacy `O/min(L_a,L_b)` overlap is not used for matching/filtering unless explicitly labeled diagnostic.
- Нет молчаливого использования sample-dependent CpG sets для статистической проверки без записи `coverage_set_mode`.
- coverage_set_sample_qc.tsv contains `total_coverage_common` and `region_sample_coverage_status`.
- coverage_set_summary.tsv records `min_region_total_coverage`.
- If `min_region_total_coverage` is used, low-coverage region × sample entries are WARN/flagged.
- There is no silent interpretation of `mu_hat_common` when `N_common == 0`.
- GLM/GLMM outputs должны включать `coverage_set_mode` и `coverage_qc_status`, когда эти поля могут быть рассчитаны.
- Нет claims о genome-wide DMR caller.
- Нет claims, что pooled delta является primary effect size.
"""
    (out_dir / "acceptance_criteria.md").write_text(text, encoding="utf-8")


def summarize_tsv(path: Path) -> str:
    if not path.exists():
        return "not written"
    try:
        df = read_table(path)
    except Exception:
        return "written, not parseable"
    if "status" in df.columns:
        counts = df["status"].astype(str).value_counts().to_dict()
        return ", ".join(f"{k}={v}" for k, v in counts.items())
    return f"{len(df)} rows"


def write_report(out_dir: Path, run_rows: list[dict], profile: str) -> None:
    sections = {
        "Overlap threshold sensitivity": summarize_tsv(out_dir / "overlap_threshold_sensitivity.tsv"),
        "Direction agreement": summarize_tsv(out_dir / "direction_agreement.tsv"),
        "Density denominator sensitivity": summarize_tsv(out_dir / "density_denominator_sensitivity.tsv"),
        "Coverage-set balance": summarize_tsv(out_dir / "coverage_set_summary.tsv"),
        "Delta weighting sensitivity": summarize_tsv(out_dir / "delta_weighting_summary.tsv"),
        "Delta bootstrap confidence intervals": summarize_tsv(out_dir / "delta_bootstrap_summary.tsv"),
        "GLM vs GLMM agreement": summarize_tsv(out_dir / "glm_glmm_status_summary.tsv"),
        "Model status audit": summarize_tsv(out_dir / "model_status_audit.tsv"),
        "Center vs interval-overlap projection": summarize_tsv(out_dir / "projection_sensitivity_summary.tsv"),
        "Random control occupancy": summarize_tsv(out_dir / "random_control_occupancy_summary.tsv"),
    }
    passed = [r for r in run_rows if r["status"] == "PASS"]
    warned = [r for r in run_rows if r["status"] == "WARN"]
    skipped = [r for r in run_rows if r["status"] == "SKIPPED"]
    failed = [r for r in run_rows if r["status"] == "FAIL"]
    lines = [
        "# Critical validation report",
        "",
        "## Profile",
        "",
        f"- profile: `{profile}`",
        f"- description: {PROFILE_DESCRIPTIONS[profile]}",
        "",
        "## What was tested",
        "",
        "This validation evaluates consistency and robustness of DMR candidates, not superiority of a caller.",
        "It checks overlap threshold sensitivity, delta direction agreement, occupancy density denominators, coverage-set balance, delta weighting sensitivity, GLM/GLMM status classes, model status flags, projection sensitivity, and random-control availability.",
        "",
        "## What passed",
        "",
        "\n".join(f"- {r['check']}" for r in passed) or "- None",
        "",
        "## What failed or was skipped",
        "",
        "\n".join(f"- {r['check']}: {r['note']}" for r in failed + skipped) or "- None",
        "",
        "## Sensitivity to overlap threshold",
        "",
        sections["Overlap threshold sensitivity"],
        "",
        "Reciprocal overlap was calculated as min(O/L_a, O/L_b). This requires sufficient overlap relative to both intervals and is stricter than O/min(L_a,L_b). Legacy min-length overlap, when present, is diagnostic only and is not used for filtering or matching.",
        "",
        "## Direction agreement",
        "",
        sections["Direction agreement"],
        "",
        "## Density denominator sensitivity",
        "",
        sections["Density denominator sensitivity"],
        "",
        "## Coverage-set balance",
        "",
        sections["Coverage-set balance"],
        "",
        "Per-sample aggregation can use different CpG sets across samples. Common-CpG filtering reduces coverage-set bias but may reduce the number of testable regions. Coverage-set QC is required before interpreting group-level delta methylation. If common-CpG mode removes many regions, GLM/GLMM validation should be interpreted cautiously.",
        "",
        "Common-CpG QC records both the number of common cytosines and region-level total common coverage. The optional parameter min_region_total_coverage can be used to require N_common_{R,s} >= N_region_min for each region × sample. In the current default lightweight validation this filter is not required unless supplied to check_coverage_set_balance.py; corresponding status columns are still reported.",
        "",
        "## Delta weighting sensitivity",
        "",
        sections["Delta weighting sensitivity"],
        "",
        "The main descriptive Delta_R is replicate-level effect size where biological replicates have equal weight. Pooled read-level delta is diagnostic only and shows sensitivity to uneven coverage distribution.",
        "",
        "## Delta bootstrap confidence intervals",
        "",
        sections["Delta bootstrap confidence intervals"],
        "",
        "Bootstrap confidence intervals were calculated for the replicate-level effect size Delta_R^{rep}. Resampling was performed at the biological-replicate level within each condition. These intervals quantify uncertainty of the descriptive effect size and do not constitute a genome-wide DMR test or FDR-controlled significance procedure. With small replicate numbers, especially 2 vs 2 designs, bootstrap intervals are diagnostic and should be interpreted cautiously.",
        "",
        "## GLM vs GLMM agreement",
        "",
        sections["GLM vs GLMM agreement"],
        "",
        "## Model status audit",
        "",
        sections["Model status audit"],
        "",
        "## Center vs interval-overlap projection",
        "",
        sections["Center vs interval-overlap projection"],
        "",
        "## Random control occupancy",
        "",
        sections["Random control occupancy"],
        "",
        "Random-control multiple testing is reported with BH-corrected bin-level empirical q-values and a global max-density profile p-value when random controls can be generated.",
        "",
        "## Direction agreement with effect-size threshold",
        "",
        "Direction agreement is evaluated only for effects exceeding configured |Delta| thresholds. Nearly zero effects are not interpreted as directional support.",
        "",
        "## Matching policy sensitivity",
        "",
        "Many-to-many matching is descriptive and can inflate pair counts. Best reciprocal or one-to-one greedy matching gives a stricter agreement estimate.",
        "",
        "## Interpretation limits",
        "",
        "- GLMM is confirmatory validation for selected regions, not genome-wide DMR calling.",
        "- q-values from DSS, methylKit and internal models are not directly equivalent.",
        "- Occupancy density is descriptive unless compared with an appropriate null/random control.",
        "- If random control is unavailable, metagene enrichment must be described cautiously.",
        "- Common-CpG mode controls only one source of bias: differences in the covered CpG/cytosine set between groups.",
        "- Delta is an effect size, not a p-value.",
        "- Pooled read-level delta is diagnostic and is not the primary effect-size definition.",
        "",
        "## Thesis-safe wording",
        "",
        "The critical validation block evaluates whether DMR/metagene conclusions are stable to overlap thresholds, projection choices, density denominators and model-status filters. It does not introduce a new DMR caller and does not prove biological function of genes.",
        "",
        "Для описательной агрегации допускается использование множества позиций I_s(R), зависящего от образца, однако такая агрегация может смешивать различие метилирования и различие набора покрытых цитозинов. Поэтому для статистической проверки выбранных регионов был добавлен coverage-set QC и режим common-CpG, в котором позиция включается в анализ только при достаточном покрытии в обеих группах. Это уменьшает риск смещения, связанного с неодинаковой структурой покрытия, но может уменьшать число тестируемых регионов.",
        "",
        "Common-CpG mode не утверждается как полное устранение всех bias; он контролирует только различие набора покрытых CpG/цитозинов между группами.",
        "",
        "Основная описательная величина Delta_R рассчитывается как replicate-level effect size, где биологические реплики имеют равный вес. Pooled read-level delta используется только как диагностическая величина, показывающая чувствительность к распределению покрытия.",
        "",
        "Bin-level empirical p-values требуют множественной коррекции или глобальной профильной статистики. Random-control occupancy является базовой null-моделью и не доказывает функциональное обогащение.",
        "",
        "Direction agreement оценивается только для эффектов, превышающих заданный порог |Delta|. Это предотвращает интерпретацию почти нулевых эффектов как направленных.",
        "",
        "Many-to-many matching используется для описательного анализа пересечений, но может завышать число пар. Best reciprocal или one-to-one greedy matching дают более строгую оценку согласованности между источниками.",
        "",
        "## Wrapper status counts",
        "",
        f"- PASS: {len(passed)}",
        f"- WARN: {len(warned)}",
        f"- SKIPPED: {len(skipped)}",
        f"- FAIL: {len(failed)}",
        "",
    ]
    (out_dir / "critical_validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    if args.list_profiles:
        for name, description in PROFILE_DESCRIPTIONS.items():
            checks = ",".join(check.name for check in checks_for_profile(name))
            print(f"{name}\t{description}\tchecks={checks}")
        return 0
    out_dir = ensure_out_dir(args.out_dir)
    run_rows: list[dict] = []
    for check in checks_for_profile(args.profile):
        check_name = check.name
        stdout = io.StringIO()
        stderr = io.StringIO()
        try:
            module = importlib.import_module(check.module)
            with redirect_stdout(stdout), redirect_stderr(stderr):
                returncode = int(module.main(["--out-dir", str(out_dir)]))
        except SystemExit as exc:
            returncode = int(exc.code or 0)
        except Exception as exc:
            returncode = 2
            stderr.write(str(exc))
        status, note = infer_status(out_dir, check_name, returncode)
        run_rows.append(
            {
                "check": check_name,
                "profile": args.profile,
                "entrypoint": check.entrypoint,
                "status": status,
                "returncode": returncode,
                "note": note,
                "stdout_tail": stdout.getvalue()[-1000:],
                "stderr_tail": stderr.getvalue()[-1000:],
            }
        )
    write_tsv(out_dir / "run_summary.tsv", run_rows)
    write_acceptance(out_dir)
    write_report(out_dir, run_rows, args.profile)
    return 1 if any(r["status"] == "FAIL" for r in run_rows) else 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
