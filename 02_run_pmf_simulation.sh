#!/bin/bash
# =============================================================================
# 02_run_pmf_simulation.sh — run umbrella sampling (min -> eql -> production) for N windows,
# at most MAX_CONCURRENT at a time. Safe to re-run: finished windows are skipped
# and interrupted ones resume.
#
#   ./02_run_pmf_simulation.sh                          # everything in pmf_config.sh
#   ./02_run_pmf_simulation.sh --prod-ns 10 --max-concurrent 6
#   ./02_run_pmf_simulation.sh --windows 6-11           # only these (or: --from 6 --to 11)
#   ./02_run_pmf_simulation.sh --dry-run                # show the plan, submit nothing
#
# Options: --windows SPEC  --from N  --to N  --prod-ns NS  --max-concurrent N
#          --min-steps N  --min2-steps N  --eql-ps PS  --eql2-ps PS   (empty = .mdp value)
#          --mode slurm|local  --dry-run  -h
# =============================================================================
set -euo pipefail
PMF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$PMF_ROOT"
source "$PMF_ROOT/pmf_config.sh"
source "$PMF_ROOT/simulation/pmf_lib.sh"
PMF_TAG="02"
# Run from PMF_Run; output goes to $PMF_DIR (5_PMF). Layout: umbN_equilibration/ = min/eql files, umbN/ = production only (prd.*, pullf-prd.xvg)

while [ $# -gt 0 ]; do
    case "$1" in
        --windows)        WINDOWS="$2"; shift 2 ;;
        --from)           WINDOW_FROM="$2"; shift 2 ;;
        --to)             WINDOW_TO="$2"; shift 2 ;;
        --prod-ns)        PROD_NS="$2"; shift 2 ;;
        --min-steps)      MIN_STEPS="$2"; shift 2 ;;
        --min2-steps)     MIN2_STEPS="$2"; shift 2 ;;
        --eql-ps)         EQL_PS="$2"; shift 2 ;;
        --eql2-ps)        EQL2_PS="$2"; shift 2 ;;
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --mode)           EXECUTION_MODE="$2"; shift 2 ;;
        --dry-run)        DRY=1; shift ;;
        -h|--help)        sed -n '2,18p' "$0"; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
done
resolve_windows conf
case "$EXECUTION_MODE" in slurm|local) ;; *) die "EXECUTION_MODE must be slurm or local" ;; esac
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || die "MAX_CONCURRENT must be a positive integer"

INPUT_DIR="$(abspath "$INPUT_DIR")"; CONF_DIR="$(abspath "$CONF_DIR")"; MDP_DIR="$(abspath "$MDP_DIR")"
PMF_MODE=prd; EXTEND_ID=""

# ── validate ─────────────────────────────────────────────────────────────────
[ -x "$GMX" ] || die "GMX not executable: $GMX"
for f in topol.top index.ndx; do [ -f "$INPUT_DIR/$f" ] || die "$INPUT_DIR/$f not found"; done
for s in min min2 eql eql2 prd; do [ -f "$MDP_DIR/$s.mdp" ] || die "$MDP_DIR/$s.mdp not found"; done
for g in ${XTC_GROUPS:-}; do grep -qE "^\[ *$g *\]" "$INPUT_DIR/index.ndx" || die "XTC_GROUPS: group '$g' not found in $INPUT_DIR/index.ndx"; done
grep -q WINDOW "$MDP_DIR/prd.mdp" || log "WARNING: no WINDOW placeholder in prd.mdp — every window would get the same reference distance"
if [ -n "$GMXLIB" ] && [ ! -d "$GMXLIB" ]; then die "GMXLIB does not exist: $GMXLIB"; fi

mapfile -t ALL < <(expand_windows "$WINDOWS")
PENDING=(); MISSING=(); DONE=0
for i in "${ALL[@]}"; do
    if [ -f "$(win_dir "$i")/.prd.done" ]; then DONE=$((DONE+1)); continue; fi
    PENDING+=("$i")
    [ -f "$(win_conf "$i")" ] || [ -f "$(win_eq_dir "$i")/start.gro" ] || MISSING+=("$(win_conf "$i")")
done
[ "${#MISSING[@]}" -eq 0 ] || die "Missing starting structures (${#MISSING[@]}), e.g. ${MISSING[0]}"
NODIST=(); for i in "${PENDING[@]}"; do [ -n "$(win_dist "$i")" ] || NODIST+=("$i"); done
[ "${#NODIST[@]}" -eq 0 ] || die "No distance in selected_configs.dat for window(s): ${NODIST[*]}"

PROD_STEPS=$(ps_to_nsteps "$(awk -v n="$PROD_NS" 'BEGIN{print n*1000}')")
log "System            : $SYS   ($EXECUTION_MODE, gmx: $GMX)"
log "Windows selected  : ${#ALL[@]}  (already finished: $DONE, to run: ${#PENDING[@]})"
log "Distance range    : $(win_dist "${ALL[0]}") -> $(win_dist "${ALL[-1]}") nm"
log "Equilibration     : min ${MIN_STEPS:-template} steps | min2 ${MIN2_STEPS:-template} steps | eql ${EQL_PS:-template} ps | eql2 ${EQL2_PS:-template} ps"
log "Production        : ${PROD_NS} ns per window  (${PROD_STEPS} steps at dt=${DT})"
if [ "$EXECUTION_MODE" = local ]; then log "Concurrency       : one umbrella at a time (local), ${NT} threads"
else log "Concurrency       : max ${MAX_CONCURRENT} umbrellas at a time, ${SLURM_NTASKS} cores each"; fi
log "Scratch           : ${SCRATCH_DIR:-none (runs directly in $PMF_DIR)}"
if [ "$EXECUTION_MODE" = slurm ]; then log "Walltime          : ${SLURM_TIME:-partition default  (WARNING: SLURM_TIME empty -> mdrun gets no -maxh)}"; fi

if [ "${#PENDING[@]}" -eq 0 ]; then log "All selected windows already finished."; exit 0; fi

mkdir -p "$PMF_OUT/logs"
SNAP="$PMF_OUT/logs/config_prd_$(date +%Y%m%d_%H%M%S).sh"
write_snapshot "$SNAP"
dispatch_windows "${PENDING[@]}"
