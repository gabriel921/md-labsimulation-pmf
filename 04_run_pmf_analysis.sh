#!/bin/bash
# =============================================================================
# 04_run_pmf_analysis.sh — WHAM over all finished windows -> 5_PMF/analysis/profile.xvg
#
#   ./04_run_pmf_analysis.sh                       # settings from pmf_config.sh
#   ./04_run_pmf_analysis.sh --begin-ps 2000       # discard first 2 ns of every window
#   ./04_run_pmf_analysis.sh --bootstrap 200       # error bars
#
# Options: --begin-ps PS  --end-ps PS  --bootstrap N  --unit kJ|kCal|kT  --temp K
# =============================================================================
set -euo pipefail
PMF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$PMF_ROOT"
source "$PMF_ROOT/pmf_config.sh"
source "$PMF_ROOT/simulation/pmf_lib.sh"
PMF_TAG="04"
WHAM_BOOTSTRAP="${WHAM_BOOTSTRAP:-0}"   # blank = no bootstrap
if [ -z "${WHAM_TEMP:-}" ]; then WHAM_TEMP="${TEMPERATURE:-}"; fi   # WHAM_TEMP empty = same as the simulations (TEMPERATURE)

while [ $# -gt 0 ]; do
    case "$1" in
        --begin-ps)  WHAM_BEGIN_PS="$2"; shift 2 ;;
        --windows)   WINDOWS="$2"; shift 2 ;;
        --from)      WINDOW_FROM="$2"; shift 2 ;;
        --to)        WINDOW_TO="$2"; shift 2 ;;
        --end-ps)    WHAM_END_PS="$2"; shift 2 ;;
        --bootstrap) WHAM_BOOTSTRAP="$2"; shift 2 ;;
        --unit)      WHAM_UNIT="$2"; shift 2 ;;
        --temp)      WHAM_TEMP="$2"; shift 2 ;;
        -h|--help)   sed -n '2,11p' "$0"; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done
[ -x "$GMX" ] || die "GMX not executable: $GMX"
[ -n "$WHAM_TEMP" ] || die "WHAM_TEMP and TEMPERATURE are both empty"
if [ -n "$GMXLIB" ]; then export GMXLIB; fi

resolve_windows existing
mapfile -t ALL < <(expand_windows "$WINDOWS")

AN="$PMF_OUT/analysis"   # results (not the plotting code in PMF_Run/analysis)
mkdir -p "$AN"; : > "$AN/tpr.dat"; : > "$AN/pullf.dat"; : > "$AN/_log_check.log"
USED=(); SKIPPED=()
printf "%-7s %-9s %-12s %s\n" window "d (nm)" "sampled(ns)" status
for i in "${ALL[@]}"; do
    d="$(win_dir "$i")"
    if compgen -G "$d/prd_ext_*.tpr" > /dev/null; then
        SKIPPED+=("$i"); printf "%-7s %-9s %-12s %s\n" "$i" "$(win_dist "$i")" "-" "extension block NOT finished — excluded (finish it with 03 --resume)"
    elif [ -f "$d/prd.tpr" ] && [ -f "$d/pullf-prd.xvg" ] && [ -f "$d/.prd.done" ]; then
        USED+=("$i")
        echo "$d/prd.tpr"        >> "$AN/tpr.dat"      # numeric order, same in both lists
        echo "$d/pullf-prd.xvg"  >> "$AN/pullf.dat"
        last=$(grep -v '^[#@]' "$d/pullf-prd.xvg" | tail -1 | awk '{printf "%.2f", $1/1000}')
        printf "%-7s %-9s %-12s ok\n" "$i" "$(win_dist "$i")" "$last"
        { echo "== window $i"; grep 'Statistics over' "$d/prd.log" 2>/dev/null || true; } >> "$AN/_log_check.log"
    else
        SKIPPED+=("$i"); printf "%-7s %-9s %-12s %s\n" "$i" "$(win_dist "$i")" "-" "NOT finished — excluded"
    fi
done
[ "${#USED[@]}" -ge 2 ] || die "Need at least 2 finished windows (found ${#USED[@]})"
[ "${#SKIPPED[@]}" -eq 0 ] || log "WARNING: ${#SKIPPED[@]} window(s) excluded: ${SKIPPED[*]}"

CMD=("$GMX" wham -it tpr.dat -if pullf.dat -temp "$WHAM_TEMP" -unit "$WHAM_UNIT" -o profile.xvg -hist histo.xvg)
if [ -n "$WHAM_BEGIN_PS" ]; then CMD+=(-b "$WHAM_BEGIN_PS"); fi
if [ -n "$WHAM_END_PS" ];   then CMD+=(-e "$WHAM_END_PS"); fi
if [ "$WHAM_BOOTSTRAP" -gt 0 ]; then CMD+=(-nBootstrap "$WHAM_BOOTSTRAP" -bsres bsResult.xvg -bsprof bsProfs.xvg); fi
if [ -n "$WHAM_ZPROF0" ];   then CMD+=(-zprof0 "$WHAM_ZPROF0"); fi
if [ -n "$WHAM_EXTRA" ];    then CMD+=($WHAM_EXTRA); fi

log "Running: ${CMD[*]}  (in $PMF_DIR/analysis/)"
( cd "$AN" && "${CMD[@]}" 2>&1 | tee wham.log )
log "PMF written to $AN/profile.xvg (histograms: $AN/histo.xvg)"
if [ "$WHAM_BOOTSTRAP" -gt 0 ]; then log "Bootstrap ($WHAM_BOOTSTRAP): average profile + standard deviations in $AN/bsResult.xvg, all profiles in $AN/bsProfs.xvg"; fi

if [ "$PLOT" = "yes" ]; then
    python3 "$PMF_ROOT/analysis/plot_pmf.py" "$AN" "$WHAM_UNIT" \
        && log "Plots: $AN/pmf.png, $AN/histograms.png" \
        || log "WARNING: plotting failed (need python3, numpy, matplotlib) — .xvg files are fine"
    if [ "$WHAM_BOOTSTRAP" -gt 0 ] && [ -f "$AN/bsResult.xvg" ]; then
        python3 "$PMF_ROOT/analysis/plot_pmf_errors.py" "$AN" "$WHAM_UNIT" \
            && log "Plot with errors: $AN/pmf_errors.png" \
            || log "WARNING: error plot failed — see $AN/bsResult.xvg"
    fi
fi
