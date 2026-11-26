#!/bin/bash
# =============================================================================
# 05_corrected_pmf.sh — Entropic correction + zero-reference shift of the WHAM
# bootstrap PMF (bsResult.xvg).
#
# Why: the PMF from gmx wham is a function of the 3D distance r between two
# centres of mass. The unbiased distribution contains the volume of the
# spherical shell, P(r) ~ r^2 exp(-G(r)/RT), so the raw profile is
#     W(r) = G(r) - 2RT ln(r)      ->      G(r) = W(r) + 2RT ln(r/r0)
# (r0 only changes a constant, which the zero-shift below removes anyway).
# Valid when the pull coordinate is a plain COM-COM "distance" (not
# "direction"/"cylinder"). Use --no-entropic to apply only the zero-shift.
#
# The whole profile is then shifted so that the mean of a "plateau" (far-apart,
# non-interacting region) is zero: the last PLATEAU_IGNORE_LAST points are
# discarded (WHAM edge effects) and the PLATEAU_AVERAGE_POINTS points before
# them are averaged.
#
# Run from anywhere (after 04_run_pmf_analysis.sh --bootstrap N):
#   ./05_corrected_pmf.sh
#   ./05_corrected_pmf.sh --no-entropic
#   ./05_corrected_pmf.sh --ignore-last 3 --average 15
#
# Input : <PMF dir>/analysis/bsResult.xvg
# Output: <PMF dir>/analysis/Corrected_PMF/bsResult_corrected.xvg
#         <PMF dir>/analysis/Corrected_PMF/pmf_corrected.png (if PLOT=yes)
#
# Settings (edit here, export as environment variables, or use the options):
#   APPLY_ENTROPIC_CORRECTION  1/0                       --entropic / --no-entropic
#   PLATEAU_IGNORE_LAST        final points excluded     --ignore-last N   (default 5)
#   PLATEAU_AVERAGE_POINTS     points averaged           --average N       (default 10)
#   PMF_REFERENCE_DISTANCE_NM  r0 in ln(r/r0), nm        --r0 NM           (default 1.0)
#   PMF_TEMP_K                 temperature (K)           --temp K  (default WHAM_TEMP/TEMPERATURE)
# The energy unit follows WHAM_UNIT and the plot switch follows PLOT
# (both from pmf_config.sh).
# =============================================================================
set -euo pipefail

PMF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PMF_ROOT"

# The shared files may use unset variables; only enforce `set -u` for this script.
set +u
source "$PMF_ROOT/pmf_config.sh"
source "$PMF_ROOT/simulation/pmf_lib.sh"
set -u

declare -F log >/dev/null || log() { echo "[$(date '+%H:%M:%S')] $*"; }
declare -F die >/dev/null || die() { echo "[ERROR] $*" >&2; exit 1; }

# ── POST-PROCESSING SETTINGS ─────────────────────────────────────────────────
APPLY_ENTROPIC_CORRECTION="${APPLY_ENTROPIC_CORRECTION:-1}"
PLATEAU_IGNORE_LAST="${PLATEAU_IGNORE_LAST:-5}"
PLATEAU_AVERAGE_POINTS="${PLATEAU_AVERAGE_POINTS:-10}"
PMF_REFERENCE_DISTANCE_NM="${PMF_REFERENCE_DISTANCE_NM:-1.0}"
INPUT_OVERRIDE=""

# ── CLI ──────────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-entropic|--no-correction) APPLY_ENTROPIC_CORRECTION=0; shift ;;
        --entropic|--correction)       APPLY_ENTROPIC_CORRECTION=1; shift ;;
        --ignore-last) PLATEAU_IGNORE_LAST="${2:?--ignore-last needs a value}"; shift 2 ;;
        --average)     PLATEAU_AVERAGE_POINTS="${2:?--average needs a value}"; shift 2 ;;
        --r0)          PMF_REFERENCE_DISTANCE_NM="${2:?--r0 needs a value}"; shift 2 ;;
        --temp)        PMF_TEMP_K="${2:?--temp needs a value}"; shift 2 ;;
        --input)       INPUT_OVERRIDE="${2:?--input needs a file}"; shift 2 ;;
        -h|--help)
            sed -n '2,/^# ====/p' "${BASH_SOURCE[0]}" | sed -n '1,/^# ====/{/^# ====/!p}' | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
done

PMF_TEMP_K="${PMF_TEMP_K:-${WHAM_TEMP:-${TEMPERATURE:-298.15}}}"
PMF_ENERGY_UNIT="${WHAM_UNIT:-kJ}"
case "$(echo "${PLOT:-yes}" | tr '[:upper:]' '[:lower:]')" in
    no|0|false|off) DO_PLOT=0 ;;
    *)              DO_PLOT=1 ;;
esac

# ── PATHS ────────────────────────────────────────────────────────────────────
# pmf_config.sh documents PMF_DIR; PMF_OUT is accepted too if your config uses it.
PMF_BASE="${PMF_OUT:-${PMF_DIR:-$PMF_ROOT/5_PMF}}"
AN="$PMF_BASE/analysis"
INPUT="${INPUT_OVERRIDE:-$AN/bsResult.xvg}"
OUT_DIR="$(dirname "$INPUT")/Corrected_PMF"
OUTPUT="$OUT_DIR/bsResult_corrected.xvg"
PLOT_FILE="$OUT_DIR/pmf_corrected.png"

[ -f "$INPUT" ] || die "Input not found: $INPUT (run 04_run_pmf_analysis.sh with --bootstrap first)"
mkdir -p "$OUT_DIR"

log "Input: $INPUT"
log "Output directory: $OUT_DIR"
log "Entropic correction enabled: $APPLY_ENTROPIC_CORRECTION"
log "Plateau: ignore last $PLATEAU_IGNORE_LAST points; average the preceding $PLATEAU_AVERAGE_POINTS points"
log "Temperature: $PMF_TEMP_K K; energy unit: $PMF_ENERGY_UNIT"

python3 - "$INPUT" "$OUTPUT" "$PLOT_FILE" \
    "$PMF_TEMP_K" "$PMF_REFERENCE_DISTANCE_NM" \
    "$PLATEAU_IGNORE_LAST" "$PLATEAU_AVERAGE_POINTS" "$PMF_ENERGY_UNIT" \
    "$APPLY_ENTROPIC_CORRECTION" "$DO_PLOT" <<'PY'
import math
import re
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
plot_path = Path(sys.argv[3])
temperature = float(sys.argv[4])
reference_distance = float(sys.argv[5])
ignore_last = int(sys.argv[6])
average_points = int(sys.argv[7])
energy_unit = sys.argv[8].strip().lower()
apply_entropic = sys.argv[9].strip().lower() in ("1", "true", "yes", "on")
do_plot = sys.argv[10] == "1"

if temperature <= 0:
    raise SystemExit("ERROR: Temperature must be greater than zero.")
if apply_entropic and reference_distance <= 0:
    raise SystemExit("ERROR: Reference distance must be greater than zero.")
if ignore_last < 0 or average_points < 1:
    raise SystemExit("ERROR: --ignore-last must be >= 0 and --average >= 1.")

# R*T in the energy unit used by gmx wham.
if energy_unit in ("kj", "kj/mol"):
    r_times_t = 0.0083144626 * temperature
    unit_label = "kJ mol^-1"
    y_unit = "kJ mol$^{-1}$"
elif energy_unit in ("kcal", "kcal/mol"):
    r_times_t = 0.0019872042 * temperature
    unit_label = "kcal mol^-1"
    y_unit = "kcal mol$^{-1}$"
elif energy_unit == "kt":
    r_times_t = 1.0          # in kT units, RT is one energy unit
    unit_label = "kT"
    y_unit = "$k_BT$"
else:
    raise SystemExit(f"ERROR: Unsupported WHAM_UNIT '{energy_unit}'. Use kJ, kCal, or kT.")

lines = input_path.read_text().splitlines(keepends=True)
records = []

for line_index, line in enumerate(lines):
    stripped = line.strip()
    if not stripped or stripped.startswith(("#", "@")):
        continue
    fields = stripped.split()
    if len(fields) < 3:
        raise SystemExit(f"ERROR: Expected at least 3 columns on line {line_index + 1}; found {len(fields)}.")
    try:
        r, pmf, unc = float(fields[0]), float(fields[1]), float(fields[2])
    except ValueError:
        raise SystemExit(f"ERROR: Non-numeric data on line {line_index + 1}.")
    if apply_entropic and r <= 0:
        raise SystemExit(f"ERROR: r must be positive for ln(r/r0); found r={r} on line {line_index + 1}.")
    records.append({"line_index": line_index, "r": r, "pmf": pmf, "unc": unc})

required = ignore_last + average_points
if len(records) < required:
    raise SystemExit(
        f"ERROR: {len(records)} data points found, but the plateau needs at least "
        f"{required} (ignore {ignore_last} + average {average_points})."
    )

# PMF_corrected(r) = PMF(r) + 2*R*T*ln(r/r0)
if apply_entropic:
    corrected = [x["pmf"] + 2.0 * r_times_t * math.log(x["r"] / reference_distance) for x in records]
    title_desc = "Entropically corrected and zero-shifted"
    curve_label = "Corrected PMF"
else:
    corrected = [x["pmf"] for x in records]
    title_desc = "Zero-shifted (no entropic correction)"
    curve_label = "Zero-shifted PMF"

# Plateau = the N points immediately before the final M points.
plateau_end = len(records) - ignore_last
plateau_start = plateau_end - average_points
plateau_mean = sum(corrected[plateau_start:plateau_end]) / average_points
final_values = [v - plateau_mean for v in corrected]

# Rebuild the XVG: only the middle (PMF) column changes.
new_value = {x["line_index"]: v for x, v in zip(records, final_values)}
title_line = re.compile(r"^\s*@\s+title\s")
out = []
for i, line in enumerate(lines):
    if i in new_value:
        f = line.split()
        f[1] = f"{new_value[i]:.8e}"
        out.append("\t".join(f) + "\n")
    elif title_line.match(line):
        out.append(f'@    title "{title_desc} bootstrap PMF"\n')
    else:
        out.append(line)
output_path.write_text("".join(out))

if do_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not found — plot skipped (python3 -m pip install matplotlib).")
        plt = None
    if plt is not None:
        r_vals = [x["r"] for x in records]
        errs = [x["unc"] for x in records]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(r_vals, final_values, yerr=errs, fmt="-", linewidth=1.5,
                    elinewidth=0.8, capsize=2, label=curve_label)
        ax.axhline(0.0, linewidth=0.9, linestyle="--")
        ax.axvspan(r_vals[plateau_start], r_vals[plateau_end - 1], alpha=0.15,
                   label=f"Zero-reference plateau ({average_points} points)")
        ax.set_xlabel("COM-COM distance (nm)")
        ax.set_ylabel(f"PMF ({y_unit})")
        ax.set_title(f"{title_desc} PMF")
        ax.legend()
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)
        print(f"Output plot: {plot_path}")

print(f"Read {len(records)} PMF data points.")
if apply_entropic:
    print(f"Entropic correction: 2*R*T*ln(r/r0), T={temperature:g} K, r0={reference_distance:g} nm.")
else:
    print("Entropic correction: DISABLED (zero-shift only).")
print(f"Plateau: rows {plateau_start + 1}-{plateau_end} "
      f"(r = {records[plateau_start]['r']:.3f}-{records[plateau_end - 1]['r']:.3f} nm; "
      f"last {ignore_last} rows ignored).")
print(f"Plateau mean before shift: {plateau_mean:.6f} {unit_label}")
print(f"Output XVG: {output_path}")
PY

log "Corrected PMF written to $OUTPUT"
