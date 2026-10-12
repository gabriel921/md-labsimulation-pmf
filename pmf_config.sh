#!/bin/bash
# =============================================================================
# pmf_config.sh — the ONLY file you normally need to edit.
# Sourced by 02_run_pmf_simulation.sh, 03_run_pmf_simulation_extend.sh and 04_run_pmf_analysis.sh.
# Lives in PMF_Run (next to the 02_/03_/04_ scripts). Relative paths are relative to THIS folder.
# Most values can be overridden for a single run from the command line
# (e.g.  ./02_run_pmf_simulation.sh --prod-ns 20 --max-concurrent 4).
# =============================================================================

# ── SYSTEM ───────────────────────────────────────────────────────────────────
SYS="pep_dimer"                                   # name used in job names
GMX="/gromacs/gromacs-2021.2/bin/gmx"             # gmx used for every umbrella stage + WHAM
GMXLIB="/home/gabrielm/Desktop/GMXLIB"           # external force-field dir (exported to every gmx call)
INPUT_DIR="1_System_Generation"                   # folder with topol.top, index.ndx (+ *.itp, *.ff)
CONF_DIR="4_Choose_cfgs"                          # starting structures, one per window (setup step 4)
CONF_PATTERN="conf%d.gro"                         # printf pattern; %d = window number (1-based)
SELECTED_FILE="4_Choose_cfgs/selected_configs.dat"  # window distances are read from here (col 5)
PMF_DIR="5_PMF"                                   # output folder: umbN/, umbN_equilibration/, logs/, analysis/ (results)
MDP_DIR="mdp"                                     # min.mdp min2.mdp eql.mdp eql2.mdp prd.mdp (use WINDOW placeholder)

# ── WINDOWS TO RUN (02_run_pmf_simulation.sh, and 04_run_pmf_analysis.sh unless --windows/--from/--to is given) ──
# Windows and their distances come from the setup (4_Choose_cfgs: conf*.gro + selected_configs.dat).
WINDOW_FROM="1"             # first window to run (e.g. 6). Empty = from the first one available
WINDOW_TO="3"               # last window to run  (e.g. 11). Empty = up to the last one available
WINDOWS=""                  # alternative: explicit list, e.g. "6-11" or "1-9,12,20-35" (overrides FROM/TO)
                            # Empty FROM/TO/WINDOWS = all windows (02: all conf*.gro in CONF_DIR; 04: all existing umbN)
WIN_PREFIX="umb"            # production folder umbN/, equilibration folder umbN_equilibration/
EQ_SUFFIX="_equilibration"

# ── SIMULATION TIMES ─────────────────────────────────────────────────────────
# A stage time left EMPTY keeps whatever nsteps is already in your .mdp template.
DT=0.002                    # ps, written into the patched .mdp whenever a time below is set
MIN_STEPS=""                # min.mdp  steepest-descent steps
MIN2_STEPS=""               # min2.mdp steepest-descent steps
EQL_PS="50"                 # eql.mdp  length (ps)
EQL2_PS="50"                # eql2.mdp length (ps)
PROD_NS=0.01                # production length PER WINDOW (ns)  -> keep it to what fits in ONE walltime (SLURM_TIME)
TEMPERATURE="310"           # K. Empty = keep the .mdp values. Set = ref-t (all groups) in eql/eql2/prd + gen-temp in eql; 04 then uses it for WHAM too
XTC_GROUPS="Protein"        # what prd.xtc saves, e.g. "Protein" or "System" (name must exist in index.ndx). Empty = .mdp (no line there = whole system)
                            # NOTE: both are fixed when a window's prd.tpr is made; changing them later does not affect windows already started or extensions

# ── UMBRELLA (pull) SETTINGS ── empty = keep what is in the .mdp files; set = applied to min, min2, eql, eql2 and prd
PULL_TYPE="umbrella"          # pull-coord1-type      e.g. umbrella
PULL_GEOMETRY="distance"      # pull-coord1-geometry  e.g. distance
PULL_K="1000.0"               # pull-coord1-k         kJ/mol/nm^2, e.g. 2000.0
PULL_RATE="0.0"               # pull-coord1-rate      0.0 for umbrella sampling
PULL_START="no"               # pull-coord1-start     no for umbrella sampling

# ── EXTENSION (03_run_pmf_simulation_extend.sh) ──────────────────────────────
# An extension = EXTEND_BLOCKS consecutive blocks of EXTEND_NS ns, per window. Every block is a separate job
# (slurm) that starts after the previous block of the SAME window finished, so each job stays short.
# Example: PROD_NS=20 in 02, then EXTEND_NS=20 and EXTEND_BLOCKS=5 here  ->  20 + 5x20 = 120 ns per window.
EXTEND_NS=0.01               # length of ONE block (ns)               (--ns)
EXTEND_BLOCKS=3             # how many blocks in a row               (--blocks)
EXTEND_WINDOWS=""           # which windows to extend, e.g. "4-6" or "1,3,8"  (--windows)   independent of WINDOWS above
EXTEND_FROM=""              # or a range: first window               (--from)
EXTEND_TO=""                #             last window                (--to)
                            # all three empty = every window that has finished production

# ── EXECUTION ────────────────────────────────────────────────────────────────
EXECUTION_MODE="slurm"      # slurm | local   (local = ONE umbrella at a time, on this machine, NT threads; MAX_CONCURRENT is ignored)
MAX_CONCURRENT=6            # slurm only: how many umbrellas (array tasks) may run at the same time
SLURM_PARTITIONS="MD32f,MD32c,QMLM,QMHM"
SLURM_NODES=1
SLURM_NTASKS=8              # cores requested per umbrella (two umbrellas share a node if it has 16 free cores)
SLURM_TIME=""               # e.g. "24:00:00". Set it: mdrun then gets -maxh and stops cleanly before the limit
SLURM_EXTRA=""              # any extra sbatch flags, e.g. "--mem=8G"
NT_MIN=8                    # mdrun -nt for the two minimisation stages
NT=8                        # mdrun -nt for eql/eql2/prd (keep <= SLURM_NTASKS)
MAXWARN_MIN2=1              # grompp -maxwarn used for min2 only (as in the original scripts)

# ── NODE-LOCAL SCRATCH ───────────────────────────────────────────────────────
# Each umbrella runs in its own folder under SCRATCH_DIR (on the node's local disk) and is copied back to
# 5_PMF/umbN and umbN_equilibration at the end and every SCRATCH_SYNC_MIN minutes (needs rsync; otherwise only at the end).
SCRATCH_DIR="/tmp"          # empty = run directly inside 5_PMF (no scratch)
SCRATCH_SYNC_MIN=0         # minutes between copy-backs while running (0 = only at the end)
MAXH_MARGIN_MIN=2          # minutes kept free at the end of the walltime for the final copy-back

# ── ANALYSIS (04_run_pmf_analysis.sh) ────────────────────────────────────────
WHAM_TEMP=310               # K 
WHAM_BEGIN_PS=0             # gmx wham -b : first time (ps) to use (discard equilibration)
WHAM_END_PS=""              # gmx wham -e : last time (ps) to use (empty = end)
WHAM_UNIT="kJ"              # kJ | kCal | kT
WHAM_BOOTSTRAP=0            # number of bootstrap rounds (0 = off)
WHAM_ZPROF0=""              # zero the PMF at this distance (nm); empty = default
WHAM_EXTRA=""               # any extra gmx wham flags
WHAM_BEGIN_PS="0"            # = --begin-ps 2000
WHAM_BOOTSTRAP=200            # = --bootstrap 200
PLOT=yes                    # yes = also write pmf.png / histograms.png (needs python3 + numpy + matplotlib)
