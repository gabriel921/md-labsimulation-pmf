#!/bin/bash
# =============================================================================
# 03_run_pmf_simulation_extend.sh — extend the production of finished windows in BLOCKS of --ns
# nanoseconds, continuing from each window's checkpoint. Slurm: every block is its own job (per window,
# starting after the previous block of that window); local: one window at a time.
#
#   ./03_run_pmf_simulation_extend.sh                          # EXTEND_* from pmf_config.sh
#   ./03_run_pmf_simulation_extend.sh --ns 10 --blocks 5       # 5 blocks of 10 ns, every finished window
#   ./03_run_pmf_simulation_extend.sh --blocks 5 --windows 4-6 # only these windows (or --from 4 --to 6)
#   ./03_run_pmf_simulation_extend.sh --resume                 # re-submit the last extension (jobs that died
#                                                              # or hit the time limit); finished blocks are skipped
# Options: --ns NS  --blocks N  --windows SPEC  --from N  --to N  --max-concurrent N
#          --mode slurm|local  --resume  --dry-run
# Blocks are numbered continuously (b1, b2, ...): every call continues after the last block used.
# =============================================================================
set -euo pipefail
PMF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$PMF_ROOT"
source "$PMF_ROOT/pmf_config.sh"
source "$PMF_ROOT/simulation/pmf_lib.sh"
PMF_TAG="03"

# The extension has its own window selection (EXTEND_*), independent of the one used by 02
WINDOWS="${EXTEND_WINDOWS:-}"; WINDOW_FROM="${EXTEND_FROM:-}"; WINDOW_TO="${EXTEND_TO:-}"
BLOCKS="${EXTEND_BLOCKS:-1}"
FIRST_BLOCK=""; EXTEND_ID=""; RESUME=0; WINDOWS_GIVEN=0
while [ $# -gt 0 ]; do
    case "$1" in
        --ns)             EXTEND_NS="$2"; shift 2 ;;
        --blocks)         BLOCKS="$2"; shift 2 ;;
        --windows)        WINDOWS="$2"; WINDOWS_GIVEN=1; shift 2 ;;
        --from)           WINDOW_FROM="$2"; shift 2 ;;
        --to)             WINDOW_TO="$2"; shift 2 ;;
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --mode)           EXECUTION_MODE="$2"; shift 2 ;;
        --resume)         RESUME=1; shift ;;
        --dry-run)        DRY=1; shift ;;
        -h|--help)        sed -n '2,15p' "$0"; exit 0 ;;
        *) die "Unknown option: $1 (see --help)" ;;
    esac
done
case "$EXECUTION_MODE" in slurm|local) ;; *) die "EXECUTION_MODE must be slurm or local" ;; esac
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || die "MAX_CONCURRENT must be a positive integer"
[[ "$BLOCKS" =~ ^[1-9][0-9]*$ ]] || die "Number of blocks must be a positive integer"

LASTFILE="$PMF_OUT/logs/last_extend"
mkdir -p "$PMF_OUT/logs"
if [ "$RESUME" = 1 ]; then
    [ -f "$LASTFILE" ] || die "No previous extension recorded ($PMF_DIR/logs/last_extend)"
    read -r FIRST_BLOCK EXTEND_NS BLOCKS LAST_WINDOWS < "$LASTFILE"
    [[ "$FIRST_BLOCK" =~ ^[0-9]+$ ]] || die "$PMF_DIR/logs/last_extend is from an older version — run 03 normally instead of --resume"
    if [ "$WINDOWS_GIVEN" = 0 ]; then WINDOWS="$LAST_WINDOWS"; WINDOW_FROM=""; WINDOW_TO=""; fi
    log "Resuming extension: $BLOCKS block(s) of ${EXTEND_NS} ns"
else
    FIRST_BLOCK="$(next_block_number)"
fi
resolve_windows existing
awk -v n="$EXTEND_NS" 'BEGIN{exit !(n>0)}' || die "Extension must be > 0 ns"

INPUT_DIR="$(abspath "$INPUT_DIR")"; CONF_DIR="$(abspath "$CONF_DIR")"; MDP_DIR="$(abspath "$MDP_DIR")"
PMF_MODE=extend
[ -x "$GMX" ] || die "GMX not executable: $GMX"

mapfile -t ALL < <(expand_windows "$WINDOWS")
PENDING=(); NOTREADY=(); DONE=0
for i in "${ALL[@]}"; do
    d="$(win_dir "$i")"
    if [ ! -f "$d/.prd.done" ] || [ ! -f "$d/prd.cpt" ]; then NOTREADY+=("$i"); continue; fi
    todo=0
    for k in $(seq 1 "$BLOCKS"); do [ -f "$d/.ext_b$((FIRST_BLOCK+k-1)).done" ] || todo=1; done
    if [ "$todo" = 1 ]; then PENDING+=("$i"); else DONE=$((DONE+1)); fi
done
[ "${#NOTREADY[@]}" -eq 0 ] || log "WARNING: windows without finished production, skipped: ${NOTREADY[*]}"
log "Windows to extend  : ${#PENDING[@]}  (already complete: $DONE)  ->  ${PENDING[*]:-none}"
log "Blocks             : b${FIRST_BLOCK}..b$((FIRST_BLOCK+BLOCKS-1)) = $BLOCKS x ${EXTEND_NS} ns = $(awk -v b="$BLOCKS" -v n="$EXTEND_NS" 'BEGIN{print b*n}') ns more per window"
if [ "$EXECUTION_MODE" = local ]; then log "Concurrency        : one window at a time (local), ${NT} threads"
else log "Concurrency        : max ${MAX_CONCURRENT} at a time per block, ${SLURM_NTASKS} cores each; each block waits for the previous block of the same window"; fi
log "Scratch            : ${SCRATCH_DIR:-none (runs directly in $PMF_DIR)}"
if [ "$EXECUTION_MODE" = slurm ]; then log "Walltime           : ${SLURM_TIME:-partition default  (WARNING: SLURM_TIME empty -> mdrun gets no -maxh)}"; fi
[ "${#PENDING[@]}" -gt 0 ] || { log "Nothing to do."; exit 0; }

SNAP="$PMF_OUT/logs/config_extend_b${FIRST_BLOCK}.sh"
write_snapshot "$SNAP"
if [ "$DRY" != 1 ]; then
    echo "$FIRST_BLOCK $EXTEND_NS $BLOCKS $WINDOWS" > "$LASTFILE"
    if [ "$RESUME" = 0 ]; then echo $((FIRST_BLOCK+BLOCKS-1)) > "$PMF_OUT/logs/extend_counter"; fi
fi
dispatch_rounds "${PENDING[@]}"
