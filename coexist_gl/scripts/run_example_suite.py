from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
from pathlib import Path
import sys

# Allow direct execution from an unpacked release without a prior editable install.
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PACKAGE_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from coexist_gl.runner import run_analysis


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _result_readme(case_id: str, summary: dict[str, object]) -> str:
    app = summary["applicability"]
    lines = [
        f"# {case_id} coexist-gl结果",
        "",
        "本文件由`scripts/run_example_suite.py`根据机器可读结果生成。案例来源、tree选择和原阳性判定见上一级`README.md`。",
        "",
        "## 观测状态",
        "",
        f"- Hosts: {app['n_hosts']}",
        f"- Single: {app['n_single']}",
        f"- Coexist (2+): {app['n_coexist']}",
        f"- Coexist fraction: {float(app['coexist_fraction']):.6g}",
        f"- Applicability: `{app['status']}`",
        "",
    ]
    if summary.get("status") != "COMPLETE":
        lines += [
            "## 结论",
            "",
            "该case未进入正式r拟合；以上状态是程序适用性结论。",
        ]
        return "\n".join(lines) + "\n"

    clustering = summary["clustering"]
    fit = summary["fit"]
    null = fit["null"]
    free = fit["free"]
    events = fit["expected_macro_edge_transitions"]
    lines += [
        "## 系统发育聚集",
        "",
        f"- Fitch steps: {clustering['fitch_steps_observed']}",
        f"- Permutations: {clustering['permutation_reps']}",
        f"- One-sided clustering p: {float(clustering['p_cluster']):.6g}",
        f"- Status: `{clustering['status']}`",
        "",
        "## Coexistence persistence拟合",
        "",
        f"- Null: r=1, A={float(null['A']):.8g}, logL={float(null['log_likelihood']):.10g}",
        f"- Free: r={float(free['r']):.8g}, A={float(free['A']):.8g}, logL={float(free['log_likelihood']):.10g}",
        f"- Delta logL: {float(fit['delta_log_likelihood']):.8g}",
        f"- Boundary LRT p: {float(fit['p_persistence_one_sided']):.6g}",
        f"- Profile-grid interval: {fit['r_profile_low']} to {fit['r_profile_high']}",
        f"- Status: `{fit['status']}`",
        "",
        "## Posterior宏观边转移期望",
        "",
        f"- S -> S: {float(events['expected_S_to_S']):.6f}",
        f"- S -> C: {float(events['expected_S_to_C_gain']):.6f}",
        f"- C -> S: {float(events['expected_C_to_S_resolution']):.6f}",
        f"- C -> C: {float(events['expected_C_to_C_persistence']):.6f}",
        "",
        "## 解释边界",
        "",
        "这是presence-conditioned、unit-edge、S/C状态模型。r是C -> C的odds multiplier；事件期望不是D/T/L reconciliation counts。原retention/AU/LOCO阳性只用于案例选择，不进入本模型。",
    ]
    return "\n".join(lines) + "\n"


def _summary_row(case_id: str, summary: dict[str, object]) -> dict[str, object]:
    app = summary["applicability"]
    row: dict[str, object] = {
        "case_id": case_id,
        "status": summary.get("status"),
        "n_hosts": app["n_hosts"],
        "n_single": app["n_single"],
        "n_coexist": app["n_coexist"],
        "coexist_fraction": app["coexist_fraction"],
        "p_cluster": "NA",
        "r_hat": "NA",
        "A_hat": "NA",
        "logL_null": "NA",
        "logL_free": "NA",
        "delta_logL": "NA",
        "p_persistence_one_sided": "NA",
        "persistence_status": "NA",
    }
    if summary.get("status") == "COMPLETE":
        fit = summary["fit"]
        row.update(
            {
                "p_cluster": summary["clustering"]["p_cluster"],
                "r_hat": fit["free"]["r"],
                "A_hat": fit["free"]["A"],
                "logL_null": fit["null"]["log_likelihood"],
                "logL_free": fit["free"]["log_likelihood"],
                "delta_logL": fit["delta_log_likelihood"],
                "p_persistence_one_sided": fit["p_persistence_one_sided"],
                "persistence_status": fit["status"],
            }
        )
    return row


def _run_one(payload: tuple[str, bool]) -> tuple[str, dict[str, object]]:
    config_text, reuse_complete = payload
    config = Path(config_text)
    case_dir = config.parent.parent
    results = case_dir / "results"
    summary_path = results / "summary.json"
    marker = results / "SUMMARY_COMPLETE.json"
    if reuse_complete and summary_path.is_file() and marker.is_file():
        summary = _read_json(summary_path)
    else:
        summary = run_analysis(config, results)
    return case_dir.name, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--examples-root", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--reuse-complete", action="store_true")
    parser.add_argument("--refresh-result-readme", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    root = args.examples_root.resolve()
    selected = set(args.case)
    configs = sorted(root.glob("*/input/family.yaml"))
    if selected:
        configs = [path for path in configs if path.parent.parent.name in selected]

    payloads = [(str(config), args.reuse_complete) for config in configs]
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            completed = list(pool.map(_run_one, payloads))
    else:
        completed = [_run_one(payload) for payload in payloads]

    rows: list[dict[str, object]] = []
    for case_id, summary in sorted(completed):
        results = root / case_id / "results"
        results.mkdir(parents=True, exist_ok=True)
        result_readme = results / "README.md"
        if args.refresh_result_readme or not result_readme.exists():
            result_readme.write_text(_result_readme(case_id, summary), encoding="utf-8")
        rows.append(_summary_row(case_id, summary))

    output = root / "EXAMPLE_SUITE_RESULTS.tsv"
    fields = list(rows[0]) if rows else ["case_id", "status"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"cases": len(rows), "summary": str(output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
