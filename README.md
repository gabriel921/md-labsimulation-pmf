# Aβ Dimer PMF Workflow (Umbrella Sampling, CHARMM36m, GROMACS)

End-to-end workflow to compute a **potential of mean force (PMF)** for the
association of two peptides with **umbrella sampling**: starts from a PDB
structure, builds the system, pulls the two peptides apart with steered MD,
picks starting structures along the pulling path, runs one umbrella
simulation per window, and finally combines everything with WHAM into a
free-energy profile.

The example system is the Aβ16–22 dimer (PDB `2MXU`, chains A and B, central
residues 17–21 = LVFFA). The pull coordinate is the distance between the
centres of mass of the central residues of each peptide.

Four scripts drive everything:

| Script | What it does |
|---|---|
| **`01_run_pmf_setup.sh`** | System generation → steered MD → extract frames → choose umbrella windows |
| **`02_run_pmf_simulation.sh`** | Runs every umbrella: minimisation → equilibration → production |
| **`03_run_pmf_simulation_extend.sh`** | Extends the production of finished umbrellas, in short blocks |
| **`04_run_pmf_analysis.sh`** | WHAM → PMF, histograms, bootstrap error bars, plots |

Steps 2–3 (umbrellas) can run as **SLURM job arrays** or **locally**, one
umbrella at a time. Both modes use the same `grompp`/`mdrun` commands.

For the background (what umbrella sampling is, why each step exists, how to
read the result) see [`TUTORIAL.md`](TUTORIAL.md).

---

## 1. Requirements

- GROMACS (tested with 2021.2)
- Bash, awk, sed
- PyMOL (only for step 01, to add ACE/NME caps to the peptides)
- CHARMM36m force field available through `GMXLIB` (see [Section 3](#3-configuration))
- SLURM (optional — only for `EXECUTION_MODE=slurm`)
- `rsync` (optional — enables periodic copy-back from node scratch)
- python3 + numpy + matplotlib (optional — only for the plots)

---

## 2. Installation

Clone the repository; no build step is needed.

```
.
├── 01_run_pmf_setup.sh
├── 02_run_pmf_simulation.sh
├── 03_run_pmf_simulation_extend.sh
├── 04_run_pmf_analysis.sh
├── pmf_config.sh                  # settings for 02, 03 and 04
├── 2MXU_16_22_dimer.pdb           # input structure
├── mdp/
│   ├── min.mdp  min2.mdp  eql.mdp  eql2.mdp  prd.mdp
├── simulation/
│   └── pmf_lib.sh                 # shared functions (never run directly)
└── analysis/
    ├── plot_pmf.py
    └── plot_pmf_errors.py
```

Run everything **from the repository root**. Step 01 creates these folders
next to the scripts:

```
1_System_Generation/   topology, box, solvent, ions, minimised system, index.ndx
2_Steered_MD/          steered MD (pull.tpr, traj_comp.xtc, ...)
3_Extract_confgs/      one .gro per SMD frame + distances.dat
4_Choose_cfgs/         chosen umbrella structures (conf1.gro, ...) + selected_configs.dat
5_PMF/                 umbrellas and results (created by 02–04)
```

---

## 3. Configuration

There are **two** places with settings:

- the **USER CONFIGURATION** block at the top of `01_run_pmf_setup.sh` (steps 1–4);
- **`pmf_config.sh`** (steps 02, 03 and 04).

Most values in `pmf_config.sh` can also be overridden for a single run from
the command line (see each step below).

### 3.1 `01_run_pmf_setup.sh`

| Variable | Meaning |
|---|---|
| `PDB` | Input PDB file (same folder as the script) |
| `SYS` | System name prefix |
| `GMX`, `PYMOL`, `GMXLIB` | Paths to `gmx`, PyMOL, and the external force-field directory |
| `START_RESI`, `END_RESI` | First/last residue of the peptide **before** capping |
| `BOX_PADDING`, `BOX_TYPE`, `ION_CONC` | Box, solvation and ionisation |
| `FF`, `WM` | Force field / water model menu indices for `pdb2gmx` (`FF=13` = CHARMM36m in the tested setup — check your `GMXLIB`) |
| `EM_STEPS` | Steepest-descent steps |
| `TEMP`, `PULL_DT`, `PULL_TIME_NS`, `PULL_RATE`, `PULL_K`, `SAVE_EVERY_PS` | Steered MD settings |
| `CENTRAL_START`, `CENTRAL_END` | Residues (original numbering) whose centre of mass defines the pull groups |
| `PULL_GRP1`, `PULL_GRP2` | Names of the index groups used for the pull coordinate |
| `N_UMBRELLAS`, `UMBRELLA_SPACING`, `FIRST_UMBRELLA_NUM` | Umbrella selection (see below) |
| `EXECUTION_MODE`, `PARTITION`, `NNODES`, `NTASKS` | `slurm` or `local`, and SLURM resources for the steered MD |
| `SMD_USE_TMP`, `SMD_TMP_BASE` | SLURM mode: run the SMD in the node's `/tmp` and copy back at the end |

**How the pulling range is set.** The peptides are pulled apart by
`PULL_RATE × PULL_TIME_NS × 1000` nm in total (e.g. 0.001 nm/ps × 2000 ps =
2 nm). Umbrellas are then picked every `UMBRELLA_SPACING` nm, so

```
number of umbrellas ≈ pulled distance / UMBRELLA_SPACING + 1
```

Make the pull long enough to give the number of umbrellas you want, and make
sure the final distance stays below half the box (the script prints a warning
if it does not — raise `BOX_PADDING` or shorten the pull). If fewer than
`N_UMBRELLAS` fit, step 4 warns you. `UMBRELLA_SPACING=0` spreads
`N_UMBRELLAS` evenly over the whole pulled range instead.

### 3.2 `pmf_config.sh`

| Block | Variables |
|---|---|
| System / folders | `SYS`, `GMX`, `GMXLIB`, `INPUT_DIR`, `CONF_DIR`, `SELECTED_FILE`, `PMF_DIR`, `MDP_DIR` |
| Windows | `WINDOW_FROM`, `WINDOW_TO`, or an explicit `WINDOWS` list like `"1-9,12,20-35"` (empty = all) |
| Times | `MIN_STEPS`, `MIN2_STEPS`, `EQL_PS`, `EQL2_PS`, `PROD_NS`, `DT` — empty keeps the value in the `.mdp` template |
| Temperature / output | `TEMPERATURE` (overrides `ref-t`/`gen-temp`), `XTC_GROUPS` (what `prd.xtc` stores) |
| Extension | `EXTEND_NS`, `EXTEND_BLOCKS`, `EXTEND_WINDOWS` or `EXTEND_FROM`/`EXTEND_TO` |
| Execution | `EXECUTION_MODE`, `MAX_CONCURRENT`, `SLURM_PARTITIONS`, `SLURM_NODES`, `SLURM_NTASKS`, `SLURM_TIME`, `SLURM_EXTRA`, `NT_MIN`, `NT` |
| Node scratch | `SCRATCH_DIR`, `SCRATCH_SYNC_MIN`, `MAXH_MARGIN_MIN` |
| Analysis | `WHAM_TEMP`, `WHAM_BEGIN_PS`, `WHAM_END_PS`, `WHAM_UNIT`, `WHAM_BOOTSTRAP`, `WHAM_ZPROF0`, `WHAM_EXTRA`, `PLOT` |

> **The values in the shipped `pmf_config.sh` are small test values**
> (`WINDOW_TO=3`, `PROD_NS=0.01`, `EQL_PS=50`, `EXTEND_NS=0.01`,
> `WHAM_BOOTSTRAP=200`). For a real calculation use all windows and
> production lengths that suit your system, e.g. `WINDOW_FROM=""`,
> `WINDOW_TO=""`, `PROD_NS=20`, `EXTEND_NS=20`, `EXTEND_BLOCKS=5`.

**The umbrella force constant and the pull settings live in the `.mdp`
files** (`mdp/*.mdp`: `pull-coord1-k`, type, geometry, ...). Edit them there.
In the `.mdp` templates the reference distance is the placeholder `WINDOW`,
which the scripts replace with each window's distance from
`selected_configs.dat`.

**Force fields and `GMXLIB`.** `GMXLIB` points GROMACS at a directory that
*contains* your force-field folders (e.g. `GMXLIB=~/my-ffs` for
`~/my-ffs/charmm36-mar2019.ff/`), so the installed GROMACS is never modified.
This repository does not ship any force fields; obtain them from an
appropriate source and respect their license and citation requirements.

---

## 4. Running the workflow

### Step 01 — Set up the system and choose the umbrellas

Edit the USER CONFIGURATION block of `01_run_pmf_setup.sh`, then:

```bash
bash 01_run_pmf_setup.sh              # everything (see EXECUTION_MODE)
```

or one stage at a time (each runs in the current shell):

```bash
bash 01_run_pmf_setup.sh system       # caps, topology, box, water, ions, EM, index groups
bash 01_run_pmf_setup.sh smd          # steered MD
bash 01_run_pmf_setup.sh extract      # one .gro per frame + distances.dat
bash 01_run_pmf_setup.sh choose       # pick umbrellas -> 4_Choose_cfgs
```

- `EXECUTION_MODE=local`: all four stages run one after the other in your shell.
- `EXECUTION_MODE=slurm`: only the steered MD is submitted (and the script
  waits for it); the other stages run in your shell. With `SMD_USE_TMP=yes`
  the SMD runs in the node's `/tmp` and is copied back at the end.

**Result:** `4_Choose_cfgs/conf1.gro … confN.gro` and
`4_Choose_cfgs/selected_configs.dat`, which lists, for every umbrella, the
original frame, its time, its real distance and the target distance. The
target distance (column 5) is the reference distance used for that window.

### Step 02 — Run the umbrellas

Edit `pmf_config.sh`, then:

```bash
./02_run_pmf_simulation.sh --dry-run                  # show the plan, submit nothing
./02_run_pmf_simulation.sh                            # run everything in pmf_config.sh
./02_run_pmf_simulation.sh --prod-ns 20 --max-concurrent 6
./02_run_pmf_simulation.sh --windows 6-11             # only some windows (or --from 6 --to 11)
```

Each window goes through `min → min2 → eql → eql2 → prd`. In SLURM mode one
**job array** is submitted (one task per umbrella), with at most
`MAX_CONCURRENT` running at once. In local mode umbrellas run one at a time
using `NT` threads.

- **Re-running is safe:** finished windows are skipped and interrupted ones
  resume from their checkpoint.
- **Time limit:** set `SLURM_TIME` (e.g. `"24:00:00"`); `mdrun` then gets
  `-maxh` and stops cleanly before the limit, leaving a checkpoint. Run the
  same command again to continue. Choose `PROD_NS` so a window fits in one
  walltime.
- **Node scratch:** with `SCRATCH_DIR=/tmp` each umbrella runs on the node's
  local disk and is copied back to `5_PMF/` at the end (and every
  `SCRATCH_SYNC_MIN` minutes if set, needs `rsync`). Set it empty to run
  directly in `5_PMF/`.

Monitor with `squeue -u $USER` and the logs in `5_PMF/logs/`.

### Step 03 — Extend the umbrellas

Production rarely converges in one go. Step 03 continues every finished
window from its checkpoint in **blocks** of `EXTEND_NS` ns. Each block is a
separate short job that starts only after the previous block of the same
window finished.

```bash
./03_run_pmf_simulation_extend.sh --ns 20 --blocks 5               # 5 x 20 ns, all finished windows
./03_run_pmf_simulation_extend.sh --blocks 5 --windows 4-6         # only these windows
./03_run_pmf_simulation_extend.sh --resume                         # redo jobs that died or hit the time limit
./03_run_pmf_simulation_extend.sh --dry-run
```

Example: `PROD_NS=20` in step 02, then `--ns 20 --blocks 5` here gives
20 + 5 × 20 = **120 ns per window**. Blocks are numbered continuously
(`b1`, `b2`, ...), so each call continues after the last block used. You can
call 03 again later to extend further.

### Step 04 — Compute the PMF

```bash
./04_run_pmf_analysis.sh                       # settings from pmf_config.sh
./04_run_pmf_analysis.sh --begin-ps 2000       # discard the first 2 ns of every window
./04_run_pmf_analysis.sh --bootstrap 200       # error bars
./04_run_pmf_analysis.sh --windows 1-30        # only some windows
```

It runs `gmx wham` over all finished windows and prints a table of which
windows were used and how many ns each sampled. Windows that are unfinished
(or have an unfinished extension block) are **excluded and listed**; finish
them first with `03 --resume`.

Outputs in `5_PMF/analysis/`:

| File | Content |
|---|---|
| `profile.xvg` | The PMF (distance vs free energy) |
| `histo.xvg` | Umbrella histograms |
| `bsResult.xvg`, `bsProfs.xvg` | Bootstrap average ± standard deviation, and all bootstrap profiles |
| `pmf.png`, `histograms.png`, `pmf_errors.png` | Plots (if `PLOT=yes`) |
| `wham.log` | `gmx wham` output |

**Always look at `histograms.png`:** neighbouring histograms must overlap,
otherwise the PMF is not reliable there (add windows or lower the spacing).

---

## 5. Output layout of `5_PMF/`

```
5_PMF/
├── umb1/                  production only: prd.mdp/tpr/cpt/gro/xtc/log, pullf-prd.xvg, pullx-prd.xvg
├── umb1_equilibration/    start.gro, min*, eql*, eql2* (everything not needed for the PMF)
├── umb2/ ...
├── logs/                  job scripts, SLURM logs, config snapshots of every run
└── analysis/              WHAM results and plots
```

---

## 6. Typical complete run

```bash
# 1. edit the USER CONFIGURATION block in 01_run_pmf_setup.sh
bash 01_run_pmf_setup.sh

# 2. edit pmf_config.sh (all windows, real PROD_NS, SLURM_TIME, ...)
./02_run_pmf_simulation.sh --dry-run
./02_run_pmf_simulation.sh

# 3. when finished, extend (repeat as needed)
./03_run_pmf_simulation_extend.sh --ns 20 --blocks 5

# 4. analysis
./04_run_pmf_analysis.sh --begin-ps 2000 --bootstrap 200
```

---

## 7. References

- Kumar, S. et al. *The weighted histogram analysis method for free-energy
  calculations on biomolecules. I. The method.* J. Comput. Chem. 13, 1011 (1992).
- Hub, J. S., de Groot, B. L., van der Spoel, D. *g_wham — A free weighted
  histogram analysis implementation including robust error and
  autocorrelation estimates.* J. Chem. Theory Comput. 6, 3713 (2010).
- GROMACS manual, section on the pull code and umbrella sampling.
