#!/usr/bin/env python3
"""Validate derived Kohler activation/deactivation process-event records."""

from __future__ import annotations

import argparse
import glob
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


INT_RE = re.compile(r"[-+]?\d+")
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
EVENT_ACTIVATION = 512
EVENT_DEACTIVATION = 1024


def ncdump(path: Path) -> str:
    try:
        result = subprocess.run(
            ["ncdump", str(path)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise RuntimeError("ncdump is required for this validator") from None
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"ncdump failed for {path}: {exc.stderr.strip()}") from exc
    return result.stdout


def read_int_variable(text: str, name: str) -> list[int]:
    return [int(value) for value in _read_variable_tokens(text, name, INT_RE)]


def read_float_variable(text: str, name: str) -> list[float]:
    return [float(value.replace("D", "E").replace("d", "e")) for value in _read_variable_tokens(text, name, FLOAT_RE)]


def _read_variable_tokens(text: str, name: str, pattern: re.Pattern[str]) -> list[str]:
    marker = f" {name} ="
    start = text.find(marker)
    if start < 0:
        marker = f"\n{name} ="
        start = text.find(marker)
    if start < 0:
        raise ValueError(f"variable {name!r} is missing")
    equals = text.find("=", start)
    end = text.find(";", equals)
    if equals < 0 or end < 0:
        raise ValueError(f"variable {name!r} has no complete data block")
    return pattern.findall(text[equals + 1 : end])


def read_mask_values(text: str, variable: str) -> list[int]:
    block_re = re.compile(rf"(?:^|\n)\s*({re.escape(variable)}(?:_\d+)?)\s*=\s*(.*?);", re.S)
    values: list[int] = []
    for _name, data in block_re.findall(text):
        values.extend(int(value) for value in INT_RE.findall(data))
    return values


def finite_positive(values: list[float], indices: list[int]) -> bool:
    return all(math.isfinite(values[idx]) and values[idx] > 0.0 for idx in indices)


def has_overlap(records: list[dict[str, int | float]], code_a: int, code_b: int) -> bool:
    keys_a = {
        (
            rec.get("time"),
            rec.get("sd_id"),
            rec.get("dm_id"),
        )
        for rec in records
        if rec["trigger_code"] == code_a
    }
    return any((rec.get("time"), rec.get("sd_id"), rec.get("dm_id")) in keys_a for rec in records if rec["trigger_code"] == code_b)


def main() -> int:
    if len(sys.argv) == 1 or sys.argv[1:] == ["--self-check"]:
        print("validate_kohler_activation_events.py: static import/self-check passed; pass --event-glob for runtime checks")
        return 0

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-glob", required=True)
    parser.add_argument("--ordinary-glob", default="SD_selected_NetCDF_*.pe*")
    parser.add_argument("--lifecycle-glob", default="SD_lifecycle_NetCDF_*")
    parser.add_argument("--require-activation", action="store_true")
    parser.add_argument("--require-deactivation", action="store_true")
    parser.add_argument("--require-significant-activation", action="store_true")
    parser.add_argument("--require-significant-deactivation", action="store_true")
    parser.add_argument("--require-kohler-context", action="store_true")
    parser.add_argument("--forbid-significant-activation", action="store_true")
    parser.add_argument("--forbid-lifecycle-activation", action="store_true")
    parser.add_argument("--check-event-mask", action="store_true")
    parser.add_argument("--check-no-phase-change-forcing", action="store_true")
    parser.add_argument("--require-condensation-overlap", action="store_true")
    parser.add_argument("--require-evaporation-overlap", action="store_true")
    args = parser.parse_args()

    event_paths = sorted(Path(p) for p in glob.glob(args.event_glob))
    if not event_paths:
        print(f"ERROR: no event files matched {args.event_glob}", file=sys.stderr)
        return 1

    records: list[dict[str, int | float]] = []
    counts: Counter[tuple[int, int]] = Counter()
    errors: list[str] = []

    for path in event_paths:
        try:
            text = ncdump(path)
            trigger_code = read_int_variable(text, "trigger_code")
            trigger_level = read_int_variable(text, "trigger_level")
            time = read_float_variable(text, "time")
            phase_pre = read_int_variable(text, "phase_state_pre")
            phase_post = read_int_variable(text, "phase_state_post")
            try:
                sd_id = read_int_variable(text, "sd_id")
                dm_id = read_int_variable(text, "dm_id")
            except ValueError:
                sd_id = [-999] * len(trigger_code)
                dm_id = [-999] * len(trigger_code)

            if not (len(trigger_code) == len(trigger_level) == len(time) == len(phase_pre) == len(phase_post)):
                errors.append(f"{path}: inconsistent core variable lengths")
                continue

            context: dict[str, list[int] | list[float]] = {}
            if args.require_kohler_context:
                for name in (
                    "kohler_rcrit_pre",
                    "kohler_rcrit_post",
                    "kohler_margin_pre",
                    "kohler_margin_post",
                ):
                    context[name] = read_float_variable(text, name)
                for name in ("activated_state_pre", "activated_state_post"):
                    context[name] = read_int_variable(text, name)

            for idx, (code, level) in enumerate(zip(trigger_code, trigger_level)):
                counts[(code, level)] += 1
                if code in (10, 11) and level == 2 and args.forbid_significant_activation:
                    errors.append(f"{path}: trigger_code={code} has forbidden trigger_level=2 at event {idx}")
                if code in (10, 11) and args.check_no_phase_change_forcing and phase_pre[idx] != phase_post[idx]:
                    errors.append(
                        f"{path}: trigger_code={code} changed phase_state {phase_pre[idx]}->{phase_post[idx]} at event {idx}"
                    )
                if code == 10 and level in (1, 2) and args.require_kohler_context:
                    if context["activated_state_pre"][idx] != 0 or context["activated_state_post"][idx] != 1:  # type: ignore[index]
                        errors.append(f"{path}: activation event {idx} does not have activated_state 0->1")
                    if context["kohler_margin_post"][idx] <= 0.0:  # type: ignore[index]
                        errors.append(f"{path}: activation event {idx} has non-positive post margin")
                if code == 11 and level in (1, 2) and args.require_kohler_context:
                    if context["activated_state_pre"][idx] != 1 or context["activated_state_post"][idx] != 0:  # type: ignore[index]
                        errors.append(f"{path}: deactivation event {idx} does not have activated_state 1->0")
                    if context["kohler_margin_pre"][idx] <= 0.0:  # type: ignore[index]
                        errors.append(f"{path}: deactivation event {idx} has non-positive pre margin")
                    if context["kohler_margin_post"][idx] > 0.0:  # type: ignore[index]
                        errors.append(f"{path}: deactivation event {idx} has positive post margin")
                records.append(
                    {
                        "time": time[idx],
                        "trigger_code": code,
                        "trigger_level": level,
                        "sd_id": sd_id[idx],
                        "dm_id": dm_id[idx],
                    }
                )

            if args.require_kohler_context:
                activation_indices = [idx for idx, code in enumerate(trigger_code) if code == 10]
                deactivation_indices = [idx for idx, code in enumerate(trigger_code) if code == 11]
                rcrit_pre = context["kohler_rcrit_pre"]  # type: ignore[assignment]
                rcrit_post = context["kohler_rcrit_post"]  # type: ignore[assignment]
                if activation_indices and not finite_positive(rcrit_pre, activation_indices):
                    errors.append(f"{path}: activation events have non-finite/non-positive kohler_rcrit_pre")
                if activation_indices and not finite_positive(rcrit_post, activation_indices):
                    errors.append(f"{path}: activation events have non-finite/non-positive kohler_rcrit_post")
                if deactivation_indices and not finite_positive(rcrit_pre, deactivation_indices):
                    errors.append(f"{path}: deactivation events have non-finite/non-positive kohler_rcrit_pre")
                if deactivation_indices and not finite_positive(rcrit_post, deactivation_indices):
                    errors.append(f"{path}: deactivation events have non-finite/non-positive kohler_rcrit_post")
        except (RuntimeError, ValueError) as exc:
            errors.append(str(exc))

    if args.require_activation and counts[(10, 1)] <= 0:
        errors.append("missing activation trigger_code=10 trigger_level=1")
    if args.require_deactivation and counts[(11, 1)] <= 0:
        errors.append("missing deactivation trigger_code=11 trigger_level=1")
    if args.require_significant_activation and counts[(10, 2)] <= 0:
        errors.append("missing activation trigger_code=10 trigger_level=2")
    if args.require_significant_deactivation and counts[(11, 2)] <= 0:
        errors.append("missing deactivation trigger_code=11 trigger_level=2")
    if args.require_condensation_overlap and not has_overlap(records, 8, 10):
        errors.append("missing same-time/same-SD condensation+activation overlap")
    if args.require_evaporation_overlap and not has_overlap(records, 9, 11):
        errors.append("missing same-time/same-SD evaporation+deactivation overlap")

    if args.check_event_mask:
        ordinary_paths = sorted(Path(p) for p in glob.glob(args.ordinary_glob) if not p.endswith(".ids"))
        if not ordinary_paths:
            errors.append(f"no ordinary/selected files matched {args.ordinary_glob}")
        else:
            bit = EVENT_ACTIVATION if (args.require_activation or args.require_significant_activation) else EVENT_DEACTIVATION
            require_sig_bit = args.require_significant_activation or args.require_significant_deactivation
            event_hits = 0
            sig_hits = 0
            phase_values: list[int] = []
            for path in ordinary_paths:
                try:
                    text = ncdump(path)
                    event_hits += sum(1 for value in read_mask_values(text, "sd_event_mask") if value & bit)
                    sig_hits += sum(1 for value in read_mask_values(text, "sd_event_sig_mask") if value & bit)
                    phase_values.extend(read_mask_values(text, "sd_phase_change_flag"))
                except RuntimeError as exc:
                    errors.append(str(exc))
            if event_hits <= 0:
                errors.append(f"sd_event_mask does not contain required bit {bit}")
            if require_sig_bit and sig_hits <= 0:
                errors.append(f"sd_event_sig_mask does not contain required significant bit {bit}")
            if not require_sig_bit and sig_hits != 0:
                errors.append(f"sd_event_sig_mask unexpectedly contains activation/deactivation bit {bit}: {sig_hits}")
            if args.check_no_phase_change_forcing and any(value not in (0, 1) for value in phase_values):
                errors.append("sd_phase_change_flag has values outside 0/1")

    if args.forbid_lifecycle_activation:
        for path_text in glob.glob(args.lifecycle_glob):
            path = Path(path_text)
            try:
                text = ncdump(path)
                if "trigger_code" in text:
                    errors.append(f"{path}: lifecycle file must not contain trigger_code")
                if "lifecycle_code" in text:
                    lifecycle_codes = read_int_variable(text, "lifecycle_code")
                    bad = [code for code in lifecycle_codes if code in (10, 11)]
                    if bad:
                        errors.append(f"{path}: lifecycle file contains activation/deactivation-like codes {bad}")
            except (RuntimeError, ValueError) as exc:
                errors.append(str(exc))

    if errors:
        print("Kohler activation event validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        print("observed_trigger_pairs=" + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(counts.items())), file=sys.stderr)
        return 1

    print(f"event_files={len(event_paths)}")
    print("trigger_pair_counts=" + ",".join(f"{code}:{level}:{count}" for (code, level), count in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
