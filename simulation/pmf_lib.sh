#!/bin/bash
# =============================================================================
# pmf_lib.sh — shared functions for the 02_/03_/04_ scripts (PMF_Run/simulation/)
# (sourced, never run directly)
#
# Folder layout inside PMF_OUT = PMF_Run/$PMF_DIR (5_PMF), per window N:
#   umbN_equilibration/   start.gro, min*, eql*, eql2*   (everything not needed for the PMF)
#   umbN/                 prd.mdp/tpr/cpt/gro/xtc/log, pullf-prd.xvg, pullx-prd.xvg
#
# With SCRATCH_DIR set, a window runs in a private folder on the node's local disk
# and is copied back to the folders above (periodically and at the end).
#
# Exit codes of a window: 0 = finished, 3 = stopped by the time limit (checkpoint saved,
# run the same command again / 03 --resume), anything else = error.
# =============================================================================
: "${DRY:=0}"
PMF_OUT="$(readlink -m -- "$PMF_ROOT/${PMF_DIR:-5_PMF}")"   # everything produced goes here
: "${WIN_PREFIX:=umb}"
: "${EQ_SUFFIX:=_equilibration}"

log()     { echo -e "[$(date '+%H:%M:%S')] ${PMF_TAG:+[$PMF_TAG] }$*"; }
die()     { echo -e "[ERROR] $*" >&2; exit 1; }
abspath() { readlink -m -- "$1"; }

win_dir()    { echo "$PMF_OUT/${WIN_PREFIX}$1"; }                   # production folder (absolute)
win_eq_dir() { echo "$PMF_OUT/${WIN_PREFIX}$1${EQ_SUFFIX}"; }       # equilibration folder (absolute)
win_conf()   { echo "$CONF_DIR/$(printf "$CONF_PATTERN" "$1")"; }    # starting structure

# Reference distance of window $1 (nm): column 5 (target) of selected_configs.dat. Empty if not found.
win_dist() {
    local f
    f="$(abspath "${SELECTED_FILE:-$CONF_DIR/selected_configs.dat}")"
    if [ -f "$f" ]; then awk -v c="$(printf "$CONF_PATTERN" "$1")" '$1==c {print $5; exit}' "$f"; fi
}

# Window numbers of the starting structures found in CONF_DIR (conf1.gro, conf2.gro, ...), one per line
list_conf_windows() {
    local dir pre suf f m
    dir="$(abspath "$CONF_DIR")"; pre="${CONF_PATTERN%%\%d*}"; suf="${CONF_PATTERN#*%d}"
    for f in "$dir"/"$pre"*"$suf"; do
        if [ -e "$f" ]; then
            m="${f##*/}"; m="${m#"$pre"}"; m="${m%"$suf"}"
            if [[ "$m" =~ ^[0-9]+$ ]]; then echo "$m"; fi
        fi
    done | sort -n
}

# Fill WINDOWS if empty: all windows of the given source (conf|existing), limited by WINDOW_FROM / WINDOW_TO
resolve_windows() {
    local all
    if [ -n "$WINDOWS" ]; then
        if [ -n "${WINDOW_FROM:-}${WINDOW_TO:-}" ]; then log "NOTE: explicit windows list given — WINDOW_FROM/--from and WINDOW_TO/--to ignored"; fi
        return 0
    fi
    if [ "$1" = conf ]; then all="$(list_conf_windows)"; else all="$(list_existing_windows)"; fi
    WINDOWS="$(echo "$all" | awk -v lo="${WINDOW_FROM:-}" -v hi="${WINDOW_TO:-}" \
        'NF && (lo=="" || $1>=lo+0) && (hi=="" || $1<=hi+0)' | paste -sd, -)"
    if [ -z "$WINDOWS" ]; then
        if [ "$1" = conf ]; then die "No starting structures ($CONF_PATTERN) in $CONF_DIR for window range ${WINDOW_FROM:-first}-${WINDOW_TO:-last}"
        else die "No umbrella folders (${WIN_PREFIX}N) for window range ${WINDOW_FROM:-first}-${WINDOW_TO:-last} in $PMF_OUT — run 02 first"; fi
    fi
}

# Window numbers that already have an umbN folder in PMF_OUT, one per line
list_existing_windows() {
    local d m
    for d in "$PMF_OUT"/"${WIN_PREFIX}"[0-9]*; do
        if [ -d "$d" ]; then
            m="${d##*/}"; m="${m#"$WIN_PREFIX"}"
            if [[ "$m" =~ ^[0-9]+$ ]]; then echo "$m"; fi
        fi
    done | sort -n
}

# "1-9,12,20-35" -> one number per line, in order, no duplicates
expand_windows() {
    echo "$1" | tr ',' '\n' | awk '
        { gsub(/[[:space:]]/, ""); if ($0 == "") next
          if ($0 ~ /^[0-9]+-[0-9]+$/) {
              split($0, a, "-"); lo = a[1]+0; hi = a[2]+0
              if (lo > hi) { print "bad window range: " $0 > "/dev/stderr"; exit 1 }
              for (k = lo; k <= hi; k++) if (!(k in s)) { s[k] = 1; print k }
          } else if ($0 ~ /^[0-9]+$/) { n = $0+0; if (!(n in s)) { s[n] = 1; print n } }
          else { print "bad window spec: " $0 > "/dev/stderr"; exit 1 } }'
}

ps_to_nsteps() { awk -v p="$1" -v dt="$DT" 'BEGIN{printf "%d", p/dt + 0.5}'; }

# ── .mdp patching ────────────────────────────────────────────────────────────
mdp_set() {   # file key value  (replace the line, or append it)
    local f="$1" k="$2" v="$3"
    if grep -qE "^[[:space:]]*${k}[[:space:]]*=" "$f"; then
        sed -i -E "s|^([[:space:]]*${k}[[:space:]]*=).*|\1 ${v}|" "$f"
    else
        echo "$k = $v" >> "$f"
    fi
}
mdp_set_temp() {  # file T — every value of ref-t (any number of tc-grps) and gen-temp
    local f="$1" t="$2"
    sed -i -E "s|^([[:space:]]*gen[-_]temp[[:space:]]*=).*|\\1 ${t}|" "$f"
    awk -v t="$t" '
        /^[[:space:]]*ref[-_]t[[:space:]]*=/ {
            line = $0; sub(/;.*/, "", line); eq = index(line, "="); key = substr(line, 1, eq-1)
            n = split(substr(line, eq+1), v, " "); if (n < 1) n = 1
            out = key "="; for (k = 1; k <= n; k++) out = out " " t; print out; next }
        { print }' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
}
make_mdp() {  # stage out dist [time_ps] [steps] — fills WINDOW; a time sets dt+nsteps, steps sets nsteps only
    local s="$1" out="$2" dist="$3" ps="${4:-}" steps="${5:-}"
    sed "s/WINDOW/$dist/g" "$MDP_DIR/$s.mdp" > "$out"
    mdp_set "$out" include "-I$INPUT_DIR"      # where grompp looks for the .itp files of topol.top (not a command-line option)
    if [ -n "${TEMPERATURE:-}" ]; then mdp_set_temp "$out" "$TEMPERATURE"; fi
    if [ "$s" = prd ] && [ -n "${XTC_GROUPS:-}" ]; then mdp_set "$out" compressed-x-grps "$XTC_GROUPS"; fi
    if [ -n "$ps" ]; then mdp_set "$out" dt "$DT"; mdp_set "$out" nsteps "$(ps_to_nsteps "$ps")"
    elif [ -n "$steps" ]; then mdp_set "$out" nsteps "$steps"; fi
}

# One equilibration stage, run inside the umbN_equilibration folder. Skipped if already done.
md_stage() {  # dir stage coords cpt maxwarn nt dist time_ps steps
    local dir="$1" s="$2" c="$3" t="$4" mw="$5" nt="$6" dist="$7" ps="$8" steps="${9:-}"
    if [ -f "$dir/.$s.done" ]; then log "  $s: done already"; return 0; fi
    log "  $s  ($(basename "$dir"))"
    (
        cd "$dir"
        make_mdp "$s" "$s.mdp" "$dist" "$ps" "$steps"
        targs=(); if [ -n "$t" ] && [ -f "$t" ]; then targs=(-t "$t"); fi
        "$GMX" grompp -f "$s.mdp" -c "$c" ${targs[@]+"${targs[@]}"} -p "$INPUT_DIR/topol.top" -n "$INPUT_DIR/index.ndx" \
            -o "$s.tpr" -po "$s-out.mdp" -maxwarn "$mw"
        "$GMX" mdrun -deffnm "$s" -nt "$nt"
        touch ".$s.done"
    )
}

# ── time limit handling ──────────────────────────────────────────────────────
walltime_hours() {   # Slurm time ("MM", "MM:SS", "HH:MM:SS", "D-HH", "D-HH:MM", "D-HH:MM:SS") -> hours
    awk -v t="$1" 'BEGIN{
        d=0; h=0; m=0; s=0
        if (t ~ /-/) { split(t, a, "-"); d=a[1]; n=split(a[2], b, ":"); h=b[1]; if(n>=2)m=b[2]; if(n>=3)s=b[3] }
        else { n=split(t, b, ":"); if(n==1){m=b[1]} else if(n==2){m=b[1]; s=b[2]} else {h=b[1]; m=b[2]; s=b[3]} }
        printf "%.4f", d*24 + h + m/60 + s/3600 }'
}
# mdrun -maxh for what is left of the job's walltime (minus the copy-back margin). Slurm mode with SLURM_TIME only.
maxh_args() {
    if [ "${EXECUTION_MODE:-}" = slurm ] && [ -n "${SLURM_TIME:-}" ]; then
        local total elapsed
        total="$(walltime_hours "$SLURM_TIME")"; elapsed=$(( $(date +%s) - ${JOB_T0:-$(date +%s)} ))
        awk -v T="$total" -v e="$elapsed" -v m="${MAXH_MARGIN_MIN:-20}" \
            'BEGIN{r = T - e/3600 - m/60; if (r < 0.05) r = 0.05; printf "-maxh %.3f", r}'
    fi
}
# Did mdrun stop early (time limit / signal) during THIS invocation? (log lines after line $2)
run_interrupted() {  # logfile lines_before
    [ -f "$1" ] && tail -n +"$(( $2 + 1 ))" "$1" | grep -qE "Run time exceeded|Received the (INT|TERM|USR1|USR2)"
}
log_lines() { if [ -f "$1" ]; then wc -l < "$1"; else echo 0; fi; }

# ── per-window work (d = production dir, e = equilibration dir; may be on scratch) ──
window_prd() {
    local i="$1" d="$2" e="$3" dist ps
    dist="$(win_dist "$i")"
    if [ -f "$d/.prd.done" ]; then log "window $i: finished already"; return 0; fi
    log "window $i  (d = $dist nm)"
    if [ ! -f "$e/start.gro" ]; then cp "$(win_conf "$i")" "$e/start.gro"; fi

    md_stage "$e" min   start.gro "" 0                "$NT_MIN" "$dist" "" "$MIN_STEPS"
    md_stage "$e" min2  min.gro   "" "$MAXWARN_MIN2"  "$NT_MIN" "$dist" "" "$MIN2_STEPS"
    md_stage "$e" eql   min2.gro  "" 0                "$NT"     "$dist" "$EQL_PS"
    md_stage "$e" eql2  eql.gro   eql.cpt 0           "$NT"     "$dist" "$EQL2_PS"

    log "  prd  (${PROD_NS} ns)"
    (
        cd "$d"
        if [ ! -f prd.tpr ]; then
            ps=$(awk -v n="$PROD_NS" 'BEGIN{print n*1000}')
            make_mdp prd prd.mdp "$dist" "$ps"
            targs=(); if [ -f "$e/eql2.cpt" ]; then targs=(-t "$e/eql2.cpt"); fi
            "$GMX" grompp -f prd.mdp -c "$e/eql2.gro" ${targs[@]+"${targs[@]}"} -p "$INPUT_DIR/topol.top" \
                -n "$INPUT_DIR/index.ndx" -o prd.tpr -po prd-out.mdp
        fi
        mh=($(maxh_args)); before="$(log_lines prd.log)"
        "$GMX" mdrun -deffnm prd -nt "$NT" -pf pullf-prd.xvg -px pullx-prd.xvg ${mh[@]+"${mh[@]}"}   # resumes from prd.cpt if present
        if run_interrupted prd.log "$before"; then
            log "  prd stopped by the time limit — checkpoint saved. Run the same command again to continue."; exit 3
        fi
        touch .prd.done
    )
}

window_extend() {   # EXTEND_ID = one block (round) of EXTEND_NS ns
    local i="$1" d="$2" ps
    if [ -f "$d/.ext_${EXTEND_ID}.done" ]; then log "window $i: block $EXTEND_ID done already"; return 0; fi
    log "window $i: extending by ${EXTEND_NS} ns (block $EXTEND_ID)"
    (
        cd "$d"
        ps=$(awk -v n="$EXTEND_NS" 'BEGIN{print n*1000}')
        # Block-specific tpr: a requeued/resumed job reuses it instead of extending a second time
        if [ ! -f "prd_ext_${EXTEND_ID}.tpr" ]; then
            "$GMX" convert-tpr -s prd.tpr -extend "$ps" -o "prd_ext_${EXTEND_ID}.tpr"
        fi
        mh=($(maxh_args)); before="$(log_lines prd.log)"
        "$GMX" mdrun -s "prd_ext_${EXTEND_ID}.tpr" -deffnm prd -cpi prd.cpt -nt "$NT" -pf pullf-prd.xvg -px pullx-prd.xvg ${mh[@]+"${mh[@]}"}
        if run_interrupted prd.log "$before"; then
            log "  stopped by the time limit — checkpoint saved. Run 03 with --resume to continue."; exit 3
        fi
        mv "prd_ext_${EXTEND_ID}.tpr" prd.tpr
        touch ".ext_${EXTEND_ID}.done"
    )
}

# ── node-local scratch ───────────────────────────────────────────────────────
have_rsync() { command -v rsync >/dev/null 2>&1; }
# Plain copies only (no -a / -p): preserving permissions, ACLs or xattrs fails between a network home and a local /tmp.
sync_dir() {   # src dst   (checkpoint first: a restart tolerates a trajectory that is longer than its checkpoint)
    if have_rsync; then
        rsync -rt --include='*.cpt' --exclude='*' "$1"/ "$2"/
        rsync -rt "$1"/ "$2"/
    else
        cp -rf "$1"/. "$2"/
    fi
}
verify_copy() {   # src dst: every file of src must exist in dst with identical content
    local f rc=0
    for f in "$1"/* "$1"/.[!.]*; do
        if [ -f "$f" ] && ! cmp -s "$f" "$2/${f##*/}"; then log "  copy-back MISMATCH: ${f##*/}"; rc=1; fi
    done
    return $rc
}
sync_back() { if [ -n "${SCR:-}" ]; then sync_dir "$SCR/d" "$FINAL_D" && sync_dir "$SCR/e" "$FINAL_E"; fi; }
scratch_finish() {   # idempotent: stop the periodic copy, copy everything back, remove the scratch folder
    local rc=0
    if [ -n "${SYNC_PID:-}" ]; then kill "$SYNC_PID" 2>/dev/null || true; SYNC_PID=""; fi
    if [ -n "${SCR:-}" ]; then
        log "  copying results back from $SCR"
        if sync_back && verify_copy "$SCR/d" "$FINAL_D" && verify_copy "$SCR/e" "$FINAL_E"; then rm -rf "$SCR"
        else log "WARNING: copy-back failed or differs — scratch kept in $SCR on $(hostname -s)"; rc=1; fi
        SCR=""
    fi
    return $rc
}

run_window() {
    local i="$1" rc u
    if [ -n "${GMXLIB:-}" ]; then export GMXLIB; fi
    JOB_T0="$(date +%s)"
    FINAL_D="$(win_dir "$i")"; FINAL_E="$(win_eq_dir "$i")"
    mkdir -p "$FINAL_D" "$FINAL_E"
    # nothing to do? (do not copy gigabytes to scratch for that)
    if [ "$PMF_MODE" = prd ] && [ -f "$FINAL_D/.prd.done" ]; then log "window $i: finished already"; return 0; fi
    if [ "$PMF_MODE" = extend ] && [ -f "$FINAL_D/.ext_${EXTEND_ID}.done" ]; then log "window $i: block $EXTEND_ID done already"; return 0; fi

    SCR=""; SYNC_PID=""; WORK_D="$FINAL_D"; WORK_E="$FINAL_E"
    if [ -n "${SCRATCH_DIR:-}" ]; then
        u="${USER:-$(id -un)}"
        SCR="$(mktemp -d "${SCRATCH_DIR%/}/${u}_${SYS}_${WIN_PREFIX}${i}_XXXXXX")" || die "cannot create a scratch folder in $SCRATCH_DIR"
        mkdir -p "$SCR/d" "$SCR/e"
        cp -rf "$FINAL_D"/. "$SCR/d"/; cp -rf "$FINAL_E"/. "$SCR/e"/
        WORK_D="$SCR/d"; WORK_E="$SCR/e"
        trap scratch_finish EXIT; trap 'exit 143' TERM INT
        log "  running in scratch $SCR on $(hostname -s)"
        if have_rsync && awk -v m="${SCRATCH_SYNC_MIN:-0}" 'BEGIN{exit !(m>0)}'; then
            ( trap - EXIT TERM INT
              while true; do sleep "$(awk -v m="$SCRATCH_SYNC_MIN" 'BEGIN{print m*60}')"; sync_back || true; done ) &
            SYNC_PID=$!
        fi
    fi

    set +e
    (
        set -e
        case "$PMF_MODE" in
            prd)    window_prd "$i" "$WORK_D" "$WORK_E" ;;
            extend) window_extend "$i" "$WORK_D" ;;
            *)      die "PMF_MODE must be prd or extend" ;;
        esac
    )
    rc=$?
    set -e
    if ! scratch_finish && [ "$rc" -eq 0 ]; then rc=1; fi
    return "$rc"
}

# ── config snapshot (the jobs source this, so command-line overrides are kept) ──
write_snapshot() {
    local f="$1" v
    {
        echo "# effective settings, $(date)"
        for v in SYS GMX GMXLIB INPUT_DIR CONF_DIR CONF_PATTERN SELECTED_FILE MDP_DIR \
                 WINDOW_FROM WINDOW_TO WINDOWS WIN_PREFIX EQ_SUFFIX \
                 DT MIN_STEPS MIN2_STEPS EQL_PS EQL2_PS PROD_NS TEMPERATURE XTC_GROUPS EXTEND_NS EXTEND_BLOCKS EXECUTION_MODE MAX_CONCURRENT \
                 SLURM_PARTITIONS SLURM_NODES SLURM_NTASKS SLURM_TIME SLURM_EXTRA \
                 SCRATCH_DIR SCRATCH_SYNC_MIN MAXH_MARGIN_MIN \
                 NT_MIN NT MAXWARN_MIN2 PMF_ROOT PMF_DIR PMF_MODE PMF_TAG EXTEND_ID; do
            printf '%s=%q\n' "$v" "${!v-}"
        done
    } > "$f"
}

# ── submit / run the windows ─────────────────────────────────────────────────
#   slurm: one job array (one task per window), at most MAX_CONCURRENT tasks at a time
#   local: one window after another, never in parallel
# Uses EXTEND_ID (the block being run, extend mode), DEP_JID (previous block's array), ROUND_SUFFIX.
dispatch_windows() {
    local wins=("$@") n=$# stamp tag i name
    stamp="$(date +%Y%m%d_%H%M%S)"; tag="${PMF_TAG:-pmf}"; name="${tag}_${EXTEND_ID:-prd}_${stamp}"
    mkdir -p "$PMF_OUT/logs"

    if [ "$DRY" = 1 ]; then
        if [ "$EXECUTION_MODE" = local ]; then log "DRY RUN — nothing run. $n window(s), local, one at a time:"
        else log "DRY RUN — nothing submitted. $n window(s), slurm, max $MAX_CONCURRENT at a time:"; fi
        for i in "${wins[@]}"; do
            printf "    window %-4s d = %-7s nm   %s/  +  %s/\n" "$i" "$(win_dist "$i")" \
                "${WIN_PREFIX}$i" "${WIN_PREFIX}$i${EQ_SUFFIX}"
        done
        return 0
    fi

    case "$EXECUTION_MODE" in
    slurm)
        local list job jid
        list="$PMF_OUT/logs/${name}.list"; printf '%s\n' "${wins[@]}" > "$list"
        job="$PMF_OUT/logs/${name}_job.sh"
        cat > "$job" << JOBEOF
#!/bin/bash
set -euo pipefail
cd "$PMF_ROOT"
source "${SNAP:?}"
source "$PMF_ROOT/simulation/pmf_lib.sh"
EXTEND_ID="$EXTEND_ID"
i=\$(sed -n "\$((SLURM_ARRAY_TASK_ID+1))p" "$list")
run_window "\$i"
JOBEOF
        local args=(--parsable --requeue -J "${SYS}_${tag}${ROUND_SUFFIX:-}" -p "$SLURM_PARTITIONS" -N "$SLURM_NODES" -n "$SLURM_NTASKS"
                    --array="0-$((n-1))%${MAX_CONCURRENT}" --chdir="$PMF_OUT"
                    -o "logs/%x_%A_%a.out" -e "logs/%x_%A_%a.err")
        if [ -n "$SLURM_TIME" ];  then args+=(-t "$SLURM_TIME"); fi
        if [ -n "${DEP_JID:-}" ]; then args+=(--dependency="aftercorr:${DEP_JID}" --kill-on-invalid-dep=yes); fi
        if [ -n "$SLURM_EXTRA" ]; then args+=($SLURM_EXTRA); fi
        jid=$(sbatch "${args[@]}" "$job"); jid="${jid%%;*}"
        LAST_JID="$jid"
        log "Submitted job array $jid${EXTEND_ID:+ (block $EXTEND_ID)}: $n window(s), max $MAX_CONCURRENT at a time${DEP_JID:+, each after its own window in job $DEP_JID} (logs: $PMF_DIR/logs/${SYS}_${tag}${ROUND_SUFFIX:-}_${jid}_*.out)"
        ;;
    local)
        local rc k=0 logf
        for i in "${wins[@]}"; do
            k=$((k+1))
            case " ${FAILED_LIST:-} " in *" $i "*) log "window $i: skipped (it failed in an earlier block)"; continue ;; esac
            logf="$PMF_OUT/logs/${name}_win${i}.log"
            log "window $i  ($k of $n, one at a time)  log: $PMF_DIR/logs/${name}_win${i}.log"
            set +e
            bash -c 'set -euo pipefail; cd "$1"; source "$2"; source "$1/simulation/pmf_lib.sh"; EXTEND_ID="$4"; run_window "$3"' \
                _ "$PMF_ROOT" "${SNAP:?}" "$i" "$EXTEND_ID" 2>&1 | tee "$logf"
            rc=${PIPESTATUS[0]}
            set -e
            if [ "$rc" -ne 0 ]; then
                FAILED_LIST="${FAILED_LIST:-} $i"
                if [ "$rc" -eq 3 ]; then log "window $i stopped by the time limit"; else log "window $i FAILED (exit $rc)"; fi
            fi
        done
        if [ -n "${FAILED_LIST:-}" ] && [ "${ROUNDS_MODE:-0}" != 1 ]; then die "window(s) not finished:${FAILED_LIST} — see $PMF_DIR/logs/"; fi
        ;;
    esac
}

# Next free block number: 1 + the highest number of extension blocks any window already has (finished or
# in progress, old or new naming) and the last number handed out (logs/extend_counter).
next_block_number() {
    local n=0 c d
    if [ -f "$PMF_OUT/logs/extend_counter" ]; then n="$(cat "$PMF_OUT/logs/extend_counter")"; fi
    for d in "$PMF_OUT/${WIN_PREFIX}"[0-9]*; do
        case "$d" in *"$EQ_SUFFIX") continue ;; esac
        if [ -d "$d" ]; then
            c=$(ls -a "$d" | grep -cE '^\.ext_.*\.done$|^prd_ext_.*\.tpr$' || true)
            if [ "$c" -gt "$n" ]; then n=$c; fi
        fi
    done
    echo $((n+1))
}

# Extension: BLOCKS consecutive blocks, ids b<FIRST_BLOCK> .. b<FIRST_BLOCK+BLOCKS-1>, for the same windows.
# slurm: every block is its own job array; task i of block k starts only after task i of block k-1 finished OK.
# local: block 1 for all windows, then block 2, ... (a window that fails is left out of later blocks).
dispatch_rounds() {
    local k
    ROUNDS_MODE=1; FAILED_LIST=""; DEP_JID=""; LAST_JID=""
    if [ "$DRY" = 1 ]; then
        log "Blocks: $BLOCKS x ${EXTEND_NS} ns  (b${FIRST_BLOCK} ... b$((FIRST_BLOCK+BLOCKS-1)))"
        EXTEND_ID="b${FIRST_BLOCK}"; dispatch_windows "$@"; return 0
    fi
    for k in $(seq 1 "$BLOCKS"); do
        EXTEND_ID="b$((FIRST_BLOCK+k-1))"; ROUND_SUFFIX="_${EXTEND_ID}"
        dispatch_windows "$@"
        if [ "$EXECUTION_MODE" = slurm ]; then DEP_JID="$LAST_JID"; fi
    done
    if [ "$EXECUTION_MODE" = local ]; then
        if [ -n "$FAILED_LIST" ]; then die "window(s) not finished:${FAILED_LIST} — run 03 with --resume"; fi
        log "All $BLOCKS block(s) finished for $# window(s)."
    fi
}
