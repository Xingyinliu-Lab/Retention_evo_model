from __future__ import annotations

import argparse
import json

from .runner import run_analysis, validate_only


def main() -> None:
    parser = argparse.ArgumentParser(prog="coexist-gl")
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--config", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--config", required=True)
    run.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.command == "validate":
        result = validate_only(args.config)
    else:
        result = run_analysis(args.config, args.output)
    print(json.dumps(result, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
