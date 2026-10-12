#!/bin/bash -e
# =============================================================================
# pmf_setup.sh
# Master setup script for PMF (umbrella sampling)
#
# Runs sequentially:
#   1) System generation  (PyMOL caps, pdb2gmx, editconf, solvate, genion, EM)
#   2) Steered MD         (grompp + mdrun with md_pull.mdp written from below)
#   3) Extract configs    (trjconv -sep + distances between the pull groups)
#   4) Choose configs     (umbrella windows picked from distances.dat -> 4_Choose_cfgs)
#
# Folder 5_PMF is created empty for the next scripts.
#
# Usage:
#   bash pmf_setup.sh                 run everything (see EXECUTION_MODE)
#   bash pmf_setup.sh system          only step 1, in this shell
#   bash pmf_setup.sh smd             only step 2, in this shell
#   bash pmf_setup.sh extract         only step 3, in this shell
#   bash pmf_setup.sh choose          only step 4, in this shell
# =============================================================================

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
PDB=2MXU_16_22_dimer.pdb                 # Input PDB file (must be in the same directory)
SYS=pep_dimer                            # System name prefix for output files
GMX=/gromacs/gromacs-2021.2/bin/gmx
PYMOL=pymol                              # PyMOL binary (ACE/NME capping)
export GMXLIB=$HOME/Desktop/GMXLIB       # GROMACS force field library path

START_RESI=16                            # First residue of the peptide (before capping)
END_RESI=22                              # Last residue of the peptide (before capping)

BOX_PADDING=1.2                          # Minimum distance from solute to box edge (nm)
BOX_TYPE=cubic                           # Box geometry
ION_CONC=0.1                             # NaCl concentration (M)
MAXWARN=1                                # Allowed grompp warnings

FF=13                                    # 13 = CHARMM36m July 2021. Check your GMXLIB
WM=1
EM_STEPS=5000                            # Steepest-descent steps

# ── STEERED MD (same values as the md_pull.mdp) ─────────────────────────────
TEMP=298                                 # ref_t for Protein and non-Protein (K)
PULL_DT=0.001                            # dt (ps)
PULL_TIME_NS=0.5                         # Length of the SMD run (ns)
PULL_RATE=0.001                          # Pulling rate (nm/ps)
PULL_K=1000                              # Force constant (kJ mol^-1 nm^-2)
SAVE_EVERY_PS=2                          # Frame saving interval, xtc (ps)

# Pull groups: the central residues of each peptide (ORIGINAL residue numbers)
CENTRAL_START=17                         # First central residue
CENTRAL_END=21                           # Last central residue
PULL_GRP1=PEP1_central                   # Index group names (also used in step 3)
PULL_GRP2=PEP2_central

# ── UMBRELLA SELECTION (step 4, from 3_Extract_confgs/distances.dat) ─────────
N_UMBRELLAS=35                           # Max number of umbrellas (e.g. 30-40)
UMBRELLA_SPACING=0.05                    # Distance between umbrellas (nm); 0 = auto: N spread evenly
FIRST_UMBRELLA_NUM=1                     # Number of the first chosen config (conf1.gro, conf2.gro, ...)

# ── EXECUTION MODE ────────────────────────────────────────────────────────────
# "slurm" — only step 2 (steered MD) is submitted to SLURM; the script waits for it to finish.
#           Steps 1, 3 and 4 run in this shell.
# "local" — steps 1, 2, 3 and 4 run one after the other in this shell.
EXECUTION_MODE=slurm
PARTITION="QMLM,MD32f"                   # SLURM partition(s)
NNODES=1                                 # Nodes per job
NTASKS=16                                # Cores per job
SMD_USE_TMP=yes                          # slurm mode: run the SMD in the node's /tmp, copy back at the end
SMD_TMP_BASE=/tmp                        # Node-local scratch directory
# ── END USER CONFIGURATION ───────────────────────────────────────────────────

# ── FOLDER STRUCTURE ─────────────────────────────────────────────────────────
SYSGEN=1_System_Generation
SMD=2_Steered_MD
EXTRACT=3_Extract_confgs
CHOOSE=4_Choose_cfgs
PMF=5_PMF

TOP=topol.top
NDX=index.ndx
ROOT=$(pwd)
SELF=$(readlink -f "$0")
PULL_PS=$(awk -v t="$PULL_TIME_NS" 'BEGIN{printf "%d", t*1000}')   # SMD length in ps
NRES=$((END_RESI - START_RESI + 1 + 2))   # residues per chain incl. ACE and NME
DIST_SEL="com of group \"$PULL_GRP1\" plus com of group \"$PULL_GRP2\""

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
log()        { echo -e "\n[$(date '+%H:%M:%S')] >>> $*"; }
die()        { echo -e "\n[ERROR] $*" >&2; exit 1; }
check_file() { [ -f "$1" ] || die "Required file not found: $1"; }
check_cmd()  { command -v "$1" &>/dev/null || die "Command not found: $1"; }
nsteps()     { awk -v t="$1" -v dt="$PULL_DT" 'BEGIN{printf "%d", t/dt}'; }   # t in ps

# Copies the SMD results from the node's /tmp back to $SMD and deletes the /tmp dir (only if the copy is OK)
smd_copy_back() {
    [ -d "$WORK" ] && [ "$WORK" != "$ROOT/$SMD" ] || return 0
    cp -df "$WORK"/* "$ROOT/$SMD"/ || { echo "[ERROR] Copy back failed — files kept in $(hostname):$WORK" >&2; return 1; }
    for f in "$WORK"/*; do
        cmp -s "$f" "$ROOT/$SMD/$(basename "$f")" || { echo "[ERROR] Copy differs for $f — files kept in $(hostname):$WORK" >&2; return 1; }
    done
    rm -rf "$WORK"
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — System generation
# ══════════════════════════════════════════════════════════════════════════════
system() {
    log "Preflight checks (step 1)..."
    check_cmd  "$GMX"
    check_cmd  "$PYMOL"
    check_file "$ROOT/$PDB"

    mkdir -p "$ROOT/$SYSGEN"; cd "$ROOT/$SYSGEN"
    cp "$ROOT/$PDB" .

    log "STEP 0 — PyMOL: adding ACE/NME caps to every chain..."
    cat > cap_pymol.py << EOF
from pymol import cmd, editor

pdb_file = "$PDB"
output_file = "peptides_capped.pdb"
start_resi = $START_RESI
end_resi   = $END_RESI

cmd.load(pdb_file, "pep")

chains = []
cmd.iterate("pep", "chains.append(chain) if chain not in chains else None", space=globals())

for chain in chains:
    n_term_selection = f"pep and chain {chain} and resi {start_resi} and name N"
    editor.attach_amino_acid(n_term_selection, "ace")

for chain in chains:
    c_term_selection = f"pep and chain {chain} and resi {end_resi} and name C"
    editor.attach_amino_acid(c_term_selection, "nme")

cmd.save(output_file, "pep")
print(f"Saved capped PDB as {output_file}")
EOF
    $PYMOL -cq cap_pymol.py
    check_file peptides_capped.pdb
    # PyMOL writes a TER before each NME cap, which makes pdb2gmx split every chain in two.
    # Keep exactly one TER at the end of each chain.
    awk '/^TER/{next} /^(ATOM|HETATM)/{c=substr($0,22,1); if(p!="" && c!=p) print "TER"; p=c} /^END/{print "TER"} {print}' \
        peptides_capped.pdb > capped_tmp.pdb && mv capped_tmp.pdb peptides_capped.pdb

    log "STEP 0 — pdb2gmx: generating topology and .gro..."
    # Capped chains need terminus type "None", otherwise pdb2gmx tries to build
    # NH3+/COO- onto ACE/NME ("atom N not found in building block 1ACE").
    # The menu number of "None" depends on the force field, so read it from pdb2gmx.
    # Probes only need the menu. pdb2gmx then waits for input forever, so each probe is killed after
    # 15 s; stdbuf makes the menu appear immediately. Output is kept in probe_*.log
    printf "%s\n%s\n" "$FF" "$WM" | timeout 15 stdbuf -o0 -e0 $GMX pdb2gmx -f peptides_capped.pdb -ignh -ter -o probe.gro 2>&1 | head -n 400 > probe_start.log || true
    TN=$(awk '/Select start terminus/{f=1;next} f && /^ *[0-9]+: None/{sub(/:.*/,""); gsub(/ /,""); print; exit}' probe_start.log)
    printf "%s\n%s\n%s\n" "$FF" "$WM" "$TN" | timeout 15 stdbuf -o0 -e0 $GMX pdb2gmx -f peptides_capped.pdb -ignh -ter -o probe.gro 2>&1 | head -n 400 > probe_end.log || true
    TC=$(awk '/Select end terminus/{f=1;next} f && /^ *[0-9]+: None/{sub(/:.*/,""); gsub(/ /,""); print; exit}' probe_end.log)
    [ -n "$TN" ] && [ -n "$TC" ] || die "Could not find the terminus type 'None' in the pdb2gmx menu — check $SYSGEN/probe_start.log and probe_end.log"
    log "Terminus 'None': start=$TN end=$TC (chains A and B)"
    printf "%s\n%s\n%s\n%s\n%s\n%s\n" "$FF" "$WM" "$TN" "$TC" "$TN" "$TC" | \
        $GMX pdb2gmx -f peptides_capped.pdb -ignh -ter -o ${SYS}.gro
    check_file "${SYS}.gro"
    check_file "$TOP"

    log "STEP 1 — editconf: building $BOX_TYPE box (padding = ${BOX_PADDING} nm)..."
    $GMX editconf \
        -f  ${SYS}.gro \
        -o  ${SYS}_box.gro \
        -c \
        -d  $BOX_PADDING \
        -bt $BOX_TYPE
    check_file "${SYS}_box.gro"

    log "STEP 2 — solvate: filling box with water..."
    $GMX solvate \
        -cp ${SYS}_box.gro \
        -p  $TOP \
        -cs spc216.gro \
        -o  ${SYS}_solv.gro
    check_file "${SYS}_solv.gro"

    cat > em.mdp << EOF
integrator      = steep
emtol           = 1000
emstep          = 0.01
nsteps          = $EM_STEPS
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

    log "STEP 3 — genion: neutralising system and adding ${ION_CONC} M NaCl..."
    $GMX grompp \
        -f em.mdp \
        -c ${SYS}_solv.gro \
        -p $TOP \
        -o ions.tpr \
        -maxwarn $MAXWARN

    echo "SOL" | $GMX genion \
        -s ions.tpr \
        -o ${SYS}_solv_ions.gro \
        -p $TOP \
        -pname SOD \
        -nname CLA \
        -neutral \
        -conc $ION_CONC
    check_file "${SYS}_solv_ions.gro"

    log "STEP 4 — energy minimisation..."
    $GMX grompp -f em.mdp -c ${SYS}_solv_ions.gro -p $TOP -o em.tpr -maxwarn $MAXWARN
    $GMX mdrun -deffnm em
    check_file em.gro

    # Index: default groups + PEP1, PEP2 (whole peptides incl. caps) and the central
    # residues of each (pull groups). Uses residue INDEX ('ri', 1-based, in file order),
    # which does not depend on the residue numbers in the .gro: ACE = 1, original
    # residue r = r - START_RESI + 2, chain B starts at NRES+1.
    log "STEP 5 — index groups..."
    C1=$((CENTRAL_START - START_RESI + 2)); C2=$((CENTRAL_END - START_RESI + 2))
    echo q | $GMX make_ndx -f em.gro -o $NDX
    n=$(grep -c '^\[' $NDX)                          # index of the next group to be created
    {
        for g in "1-$NRES PEP1" "$((NRES+1))-$((2*NRES)) PEP2" \
                 "$C1-$C2 $PULL_GRP1" "$((C1+NRES))-$((C2+NRES)) $PULL_GRP2"; do
            set -- $g
            echo "ri $1"; echo "name $n $2"; n=$((n+1))
        done
        echo q
    } | $GMX make_ndx -f em.gro -n $NDX -o ${NDX}_new.ndx
    mv ${NDX}_new.ndx $NDX
    for g in PEP1 PEP2 $PULL_GRP1 $PULL_GRP2; do
        grep -q "^\[ $g \]" $NDX || die "Index group $g not found in $SYSGEN/$NDX"
    done
    awk '/^\[/{if(g!="")print "  " g ": " n " atoms"; g=$2; n=0; next} {n+=NF} END{print "  " g ": " n " atoms"}' $NDX | \
        grep -E "PEP1|PEP2"
    cd "$ROOT"
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Steered MD
# ══════════════════════════════════════════════════════════════════════════════
smd() {
    log "Preflight checks (step 2)..."
    check_cmd  "$GMX"
    for f in em.gro $TOP $NDX; do check_file "$ROOT/$SYSGEN/$f"; done

    mkdir -p "$ROOT/$SMD"; cd "$ROOT/$SMD"
    cp "$ROOT/$SYSGEN"/{em.gro,$TOP,$NDX} .
    cp "$ROOT/$SYSGEN"/*.itp . 2>/dev/null || true

    cat > md_pull.mdp << EOF
title       = Harmonic pulling simulation
;define      = -DPOSRES

; Run parameters
integrator  = md
dt          = $PULL_DT
tinit       = 0
nsteps      = $(nsteps $PULL_PS)

; COM motion removal
nstcomm         = 10
comm-mode       = Linear
comm-grps       = System

; Energy groups
energygrps = PEP1 PEP2 SOL

; Output parameters
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstenergy   = 1000
nstlog      = 1000

nstxout-compressed = $(nsteps $SAVE_EVERY_PS)
compressed-x-precision = 1000

; Bond parameters
constraint_algorithm    = lincs
constraints             = all-bonds
lincs_iter              = 2
lincs_order             = 4
continuation            = yes

; Neighborsearching
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
verlet-buffer-tolerance = 1.0e-06

; PME electrostatics parameters
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
ewald_rtol      = 1e-5

; Temperature coupling
tcoupl      = V-rescale
tc_grps     = Protein  non-Protein
tau_t       = 0.1  0.1
ref_t       = $TEMP  $TEMP

; Generate velocities is off
gen_vel     = no

; Periodic boundary conditions are on in all directions
pbc     = xyz

; Long-range dispersion correction
DispCorr    = EnerPres

; Pull code
pull                    = yes
pull_ngroups            = 2
pull_ncoords            = 1
pull_group1_name        = $PULL_GRP1
pull_group2_name        = $PULL_GRP2
pull_coord1_type        = umbrella
pull_coord1_geometry    = distance
pull_coord1_groups      = 1 2
pull_coord1_dim         = Y Y Y
pull_coord1_rate        = $PULL_RATE
pull_coord1_k           = $PULL_K
pull_coord1_start       = yes
EOF

    # The pulled distance must stay below half the box
    d0=$($GMX distance -s em.gro -f em.gro -n $NDX -select "$DIST_SEL" -oall d0.xvg >/dev/null 2>&1 \
         && awk '!/^[#@]/{print $2; exit}' d0.xvg || true)
    awk -v d0="${d0:-0}" -v r=$PULL_RATE -v t=$PULL_PS 'END{ \
        if (d0 + r*t > $1/2) print "[WARNING] Final pull distance " d0 + r*t " nm exceeds half the box (" $1/2 " nm): raise BOX_PADDING or lower PULL_TIME_NS/PULL_RATE" }' em.gro

    # Slurm mode: run in the node's /tmp (not on /home) and copy everything back at the end.
    # A requeued job restores state.cpt/pull.tpr from $SMD (written there by the copy back below).
    WORK="$ROOT/$SMD"
    if [ "$EXECUTION_MODE" = slurm ] && [ "$SMD_USE_TMP" = yes ] && [ -n "$SLURM_JOB_ID" ]; then
        WORK="$SMD_TMP_BASE/${USER}_PMF_SMD_$SLURM_JOB_ID"
        log "SMD runs in $(hostname):$WORK (copied back to $SMD at the end)"
        mkdir -p "$WORK"
        cp md_pull.mdp em.gro $TOP $NDX "$WORK"/
        cp *.itp "$WORK"/ 2>/dev/null || true
        for f in state.cpt pull.tpr; do if [ -f "$f" ]; then cp "$f" "$WORK"/; fi; done
        trap 'exit 143' TERM                        # on scancel/requeue: run the EXIT trap below
        trap smd_copy_back EXIT                     # also copies back if mdrun fails
    fi
    cd "$WORK"

    log "STEP 6 — steered MD..."
    if [ -f state.cpt ]; then                       # requeued job: continue
        $GMX mdrun -s pull.tpr -cpi state.cpt
    else
        $GMX grompp -f md_pull.mdp -c em.gro -p $TOP -n $NDX -o pull.tpr -maxwarn $MAXWARN
        $GMX mdrun -s pull.tpr
    fi
    smd_copy_back; trap - EXIT TERM
    cd "$ROOT"
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Extract configurations
# ══════════════════════════════════════════════════════════════════════════════
extract() {
    log "Preflight checks (step 3)..."
    check_cmd "$GMX"
    for f in pull.tpr traj_comp.xtc $NDX; do check_file "$ROOT/$SMD/$f"; done

    mkdir -p "$ROOT/$EXTRACT"; cd "$ROOT/$EXTRACT"
    ln -sf "$ROOT/$SMD"/{pull.tpr,traj_comp.xtc,$NDX} .

    log "STEP 7 — trjconv: one .gro per frame (conf0.gro, conf1.gro, ...)..."
    echo -e "0\n" | $GMX trjconv -s pull.tpr -f traj_comp.xtc -o conf.gro -sep

    log "STEP 8 — distances between $PULL_GRP1 and $PULL_GRP2..."
    $GMX distance -s pull.tpr -f traj_comp.xtc -n $NDX -select "$DIST_SEL" -oall dist.xvg
    { echo "# conf time(ps) distance(nm)"
      awk '!/^[#@]/{printf "%d %s %s\n", n++, $1, $2}' dist.xvg; } > distances.dat
    cd "$ROOT"
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Choose configurations (umbrella windows)
# ══════════════════════════════════════════════════════════════════════════════
choose() {
    log "Preflight checks (step 4)..."
    check_file "$ROOT/$EXTRACT/distances.dat"

    mkdir -p "$ROOT/$CHOOSE"; cd "$ROOT/$CHOOSE"
    rm -f conf*.gro selected_configs.dat            # a re-run never mixes old and new umbrellas

    log "STEP 9 — choosing up to $N_UMBRELLAS umbrellas (spacing $UMBRELLA_SPACING nm)..."
    # First frame = first umbrella. Then, for each target distance (d_first + k*spacing), the first frame
    # that reaches it (or the one just before it, if closer) after the previously chosen frame.
    SEL=$(awk -v n=$N_UMBRELLAS -v sp=$UMBRELLA_SPACING '
        !/^#/ { m++; conf[m]=$1; time[m]=$2; dist[m]=$3; if (m==1 || $3>dmax) dmax=$3 }
        END {
            d0 = dist[1]
            if (sp <= 0) sp = (dmax - d0)/n                           # auto spacing (last target stays reachable)
            printf "# range: %.3f-%.3f nm   spacing used: %.4f nm\n", d0, dmax, sp
            k = 1; sel[1] = 1; tgt[1] = d0; prev = 1
            for (j = 2; j <= n; j++) {
                target = d0 + (j-1)*sp
                if (target > dmax + 1e-9) break
                c = 0
                for (i = prev+1; i <= m; i++) if (dist[i] >= target - 1e-9) { c = i; break }
                if (!c) break
                pick = c
                if (c-1 > prev && target - dist[c-1] < dist[c] - target) pick = c-1
                sel[++k] = pick; tgt[k] = target; prev = pick
            }
            for (s = 1; s <= k; s++) printf "%d %s %s %.3f\n", conf[sel[s]], time[sel[s]], dist[sel[s]], tgt[s]
        }' "$ROOT/$EXTRACT/distances.dat")

    {
        echo "# Umbrellas chosen from $EXTRACT/distances.dat"
        echo "$SEL" | grep '^#'
        echo "# new_conf  orig_conf  time(ps)  distance(nm)  target(nm)"
        num=$FIRST_UMBRELLA_NUM
        while read -r orig t d tgt; do
            check_file "$ROOT/$EXTRACT/conf${orig}.gro"
            cp "$ROOT/$EXTRACT/conf${orig}.gro" conf${num}.gro
            printf "conf%-6s conf%-6s %9s %9s %9s\n" "${num}.gro" "${orig}.gro" "$t" "$d" "$tgt"
            num=$((num+1))
        done <<< "$(echo "$SEL" | grep -v '^#')"
    } > selected_configs.dat

    found=$(grep -vc '^#' selected_configs.dat)
    log "$found umbrellas copied to $CHOOSE/ (list in $CHOOSE/selected_configs.dat)"
    [ "$found" -ge "$N_UMBRELLAS" ] || echo "[WARNING] Only $found of $N_UMBRELLAS umbrellas fit in the pulled range with this spacing: lower UMBRELLA_SPACING (0 = auto), or pull further (PULL_TIME_NS/PULL_RATE)" >&2
    cd "$ROOT"
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
mkdir -p $SYSGEN $SMD $EXTRACT $CHOOSE $PMF
echo "  $SYSGEN/  — system generation (caps, topology, box, solvation, ions, EM)"
echo "  $SMD/  — steered MD"
echo "  $EXTRACT/  — extracted configurations + distances.dat"
echo "  $CHOOSE/  — chosen umbrella configs + selected_configs.dat"
echo "  $PMF/  — (next script)"

case "${1:-all}" in
    system)  system ;;
    smd)     smd ;;
    extract) extract ;;
    choose)  choose ;;
    all)
        case "$EXECUTION_MODE" in
            slurm)
                system
                log "Submitting step 2 (steered MD) to SLURM and waiting for it to finish..."
                sbatch --wait --requeue -p $PARTITION -N $NNODES -n $NTASKS --chdir=$ROOT \
                    -J pmf_2_smd -o SMD.sout -e SMD.serr --wrap="bash $SELF smd"
                extract; choose ;;
            local)
                system; smd; extract; choose ;;
            *)  die "Unknown EXECUTION_MODE='$EXECUTION_MODE' — must be 'slurm' or 'local'." ;;
        esac ;;
    *)  die "Usage: bash $0 [system|smd|extract|choose]" ;;
esac
