# A Student's Guide to This Umbrella Sampling / PMF Workflow

This document explains *why* each step of this repository exists — the
physics behind it — and *how to read* the result. It assumes you know basic
molecular biology and a little physical chemistry, but nothing about free-energy
calculations. If you have not seen the all-atom workflow of this series, you
do not need to: everything needed is explained here.

---

## 1. What question are we asking?

Two peptides floating in water either stay apart or stick together. The
interesting questions are quantitative:

- **How strongly** do they stick? (How much free energy do you gain by binding?)
- **At what distance** do they bind best?
- Is there a **barrier** on the way in or out?

Our example is a dimer of the Aβ16–22 fragment (sequence KLVFFAE, from the
Alzheimer's-related Aβ peptide, taken from the fibril structure `2MXU`). Its
central residues 17–21 (**LVFFA**) are the hydrophobic core that drives Aβ
self-association — the same region that matters in the full-length peptide.
Two such fragments coming together is the simplest model of the first step of
aggregation: **dimerisation**.

The answer to all three questions is one curve: the **potential of mean
force (PMF)**.

---

## 2. What is a PMF?

Pick one number that describes the process — here, the **distance `d`**
between the two peptides (more precisely, between the centres of mass of
their central residues). Now ask: *how favourable is it for the system to be
at distance `d`, after averaging over everything else* (all the water
molecules, all the side-chain wiggling, all the orientations)?

That averaged free energy as a function of `d` is the PMF, `W(d)`. It is
directly linked to probability:

```
W(d) = − kT · ln P(d)        (up to an arbitrary constant)
```

where `P(d)` is how often the system would be found at distance `d`. High
probability = low free energy. So:

- a **well (minimum)** in the PMF = a bound state;
- the **depth** of the well relative to the far-away plateau ≈ binding free energy;
- a **hump** between them = a free-energy barrier.

Reading it is like reading a hiking map: the PMF is the altitude profile of
the path from "far apart" to "bound".

---

## 3. Why can't we just run a long simulation?

In principle, `P(d)` could be measured by running MD and counting how often
each distance occurs. In practice this fails for exactly the part we care
about. If the bound state is, say, 30 kJ/mol more stable than the unbound
one, the system will sit in the bound state and almost never visit the
barrier region or the far-apart region. You would wait a very long time to
sample the rare distances — and the free energy there is determined by those
rare visits. This is the **sampling problem** (a *rare-event* problem).

**Umbrella sampling** is a way around it: instead of waiting for the system to
visit every distance, we *force* it to stay near each distance in turn, and
then mathematically remove the effect of our forcing.

---

## 4. Umbrella sampling in plain words

Imagine mapping the altitude along a mountain path in thick fog. You cannot
walk the whole path in one go, because every time you try you slide into the
nearest valley and stay there. So instead you do this:

1. Tie yourself with a **spring (rope)** to a pole planted at kilometre 0.
   You can wander a little around the pole but not far. You record the
   altitude profile of that small stretch.
2. Move the pole to kilometre 0.05. Wander around it. Record.
3. Repeat, moving the pole in small steps along the whole path.
4. Finally, combine all the small, overlapping stretches into one profile,
   correcting for the fact that the rope was pulling you.

That is umbrella sampling:

| Mountain | Simulation |
|---|---|
| Position along the path | Distance `d` between the peptides |
| The pole | The **reference distance** `d₀` of a window |
| The spring | A **harmonic bias potential** `½ k (d − d₀)²` |
| One stretch around one pole | One **umbrella window** (an independent simulation) |
| Combining the stretches | **WHAM** |

Each window is a normal MD simulation with one extra term: a spring that
keeps `d` close to that window's `d₀`. The name "umbrella" comes from the
fact that the bias holds the system under a "shelter" around `d₀`.

Because the spring is known exactly, its effect can be subtracted later. What
is left is the unbiased free energy.

---

## 5. The workflow, step by step

### Overview

```
01  system → steered MD → extract frames → choose windows
02  run every umbrella (min → eql → eql2 → production)
03  (optional, repeatable) extend the umbrellas
04  WHAM → PMF + error bars + plots
```

### Step 01 — Building the system

**Capping (PyMOL).** The 7-residue fragment was cut out of a long peptide. A
real peptide would be continuous at both ends, but a bare fragment would end
in a charged NH₃⁺ and COO⁻ that do not exist in the full-length molecule and
would artificially attract each other. We therefore add neutral **ACE** (N-terminal
acetyl) and **NME** (C-terminal N-methylamide) caps, which mimic the missing
neighbouring residues. At `pdb2gmx` the terminus type "None" is then chosen
so GROMACS does not try to add charged ends on top of the caps (the script
finds this menu entry automatically).

**Topology, box, water, ions.** As in any MD setup: assign the **CHARMM36m**
force field (parameterised to behave well for disordered and flexible
peptides), put the dimer in a cubic box with generous padding, fill it with
explicit water, and add NaCl to neutralise the system and reach 0.1 M. The box
needs to be large — see the pulling step: the peptides are going to be pulled
apart inside it, and the final separation must stay below **half the box
length**, otherwise a peptide would start to interact with the periodic image
of the other one.

**Energy minimisation.** Removes bad contacts created by solvation.

**Index groups.** The script creates the named groups that the rest of the
workflow uses: the two whole peptides (`PEP1`, `PEP2`) and the **central
residues** of each (`PEP1_central`, `PEP2_central`). The distance between
the centres of mass of these two central groups is our **reaction
coordinate** — the single number whose free energy profile we want.

> **Choosing a reaction coordinate** is the most important scientific
> decision in umbrella sampling. A centre-of-mass distance between the
> hydrophobic cores is a natural, simple choice for dimerisation, but it
> ignores other things that might matter (relative orientation, which
> residues touch, etc.). See Section 7.

### Step 01 (continued) — Steered MD: making the starting structures

Umbrella sampling needs a **starting structure for every window**, with the
peptides already at roughly the right distance. We get them by *pulling*.

In steered MD (SMD) a spring is attached to the two groups and its reference
point moves at a constant, slow rate (here 0.001 nm/ps). The spring drags the
two peptides apart, and we save frames along the way. The result is a movie of
the dimer being pulled from bound to separated.

This pulling is a **non-equilibrium** process: we dragged the system and it
did not have time to relax. The structures along the way are therefore *not*
reliable on their own, and we must **never** use the pulling curve itself as
a free energy. That is why the SMD is only used as a way of generating
starting points. The windows then do the real work: each frame is
equilibrated at its own distance, and all the physics that matters is measured
there.

**Extract and choose.** Every saved frame is written as a `.gro` file and the
distance between the groups is computed (`distances.dat`). The script then
selects the first frame and one frame every `UMBRELLA_SPACING` nm (0.05 nm in
our example) up to `N_UMBRELLAS`, renumbers them `conf1.gro`, `conf2.gro`, ...
and records the choice in `selected_configs.dat`. The **target distance**
listed there is the `d₀` that the umbrella spring will be centred on, which is
a round, evenly spaced number rather than whatever the frame happened to
have.

*How close should the windows be?* Close enough that neighbouring windows
**overlap**, as explained below. 0.05 nm spacing with a force constant of a
couple of thousand kJ mol⁻¹ nm⁻² is a typical, safe combination.

### Step 02 — Running the umbrellas

Every window is an independent simulation, which is why they can run in
parallel on a cluster (one SLURM array task per window). Each goes through:

1. **`min`, `min2`** — two energy minimisations with the umbrella spring active.
2. **`eql`** — short NVT equilibration (temperature coupled, random initial
   velocities) with the spring on.
3. **`eql2`** — NPT equilibration (pressure also coupled) so the box finds the
   right density.
4. **`prd`** — production: this is the part that is analysed.

All the minimisation and equilibration files go into `umbN_equilibration/`,
and only production files stay in `umbN/`, which keeps the analysis folder
clean.

**Why equilibrate with the spring on?** Because the spring is part of the
physics of that window. We want the system relaxed *under the bias*, not
relaxed freely and then suddenly held.

### Step 03 — Extending the umbrellas

Free-energy results converge slowly. Instead of guessing the right simulation
length in advance, we run a first production, look at the PMF, and then
**extend**: each window continues from its checkpoint, in short blocks, so a
long simulation can be spread across many short, queue-friendly jobs.
Scientifically, ten 20 ns blocks are equivalent to one 200 ns run.

### Step 04 — WHAM: putting the windows together

After production, for each window we have a record of the distance over time
(`pullx`) and the force of the spring (`pullf`). The distances in one window
form a **histogram**, narrow and centred near that window's `d₀` (the spring
stops the system from wandering far).

The **Weighted Histogram Analysis Method (WHAM)** does two things at once:

1. It **removes the bias** of each spring from its histogram (we know the
   spring exactly, so we can undo it).
2. It **stitches** the histograms of all windows into a single, consistent
   `P(d)` and therefore a single `W(d)`. Where two windows overlap, the
   method uses that overlap to fix the free-energy offset between them.

The output is the PMF. For error bars, **bootstrapping** is used: the
analysis is repeated many times on randomly resampled data, and the spread of
the resulting curves gives a standard deviation at each distance.

---

## 6. Checking that the result is trustworthy

Getting a smooth curve out of WHAM is easy. Getting a **correct** one takes
checks. Do these every time.

**1. Histogram overlap (`histograms.png`).** Every neighbouring pair of
histograms should clearly overlap. If there is a gap, WHAM has no information
to link the windows and the PMF is unreliable there. Fix: more windows or
smaller spacing, or a softer spring (smaller `k`) so that histograms are wider.

**2. Discard the equilibration (`--begin-ps`).** The first part of each
production run still remembers the starting structure. Discarding it (e.g.
the first 2 ns) removes this bias. Rerun the analysis with different cut-offs
and see whether the PMF changes.

**3. Convergence with simulation time.** Compute the PMF using only the first
half of the data, then all of it (`--end-ps`). If the curve still changes,
**extend** (step 03) and repeat. A converged PMF stops changing when you add
more sampling.

**4. Error bars.** Use the bootstrap. Be cautious: bootstrap errors assume
the data are independent samples, so they tend to be *optimistic* for
strongly correlated data. They are a minimum estimate of uncertainty, not a
guarantee.

**5. Starting-structure dependence.** For flexible peptides, windows can get
"stuck" in different conformations depending on how they were started. If you
can, check that windows built from a *different* pulling run (or a different
random seed) give a similar PMF.

---

## 7. Reading the PMF and knowing its limits

A typical dimerisation PMF has:

- a **plateau at large distance** — the peptides no longer feel each other.
  The free energy there is conventionally set as the reference (zero). Note that for a
  1D distance coordinate, the far region also contains an *entropic* term that
  grows with distance (more volume is available), so be careful when comparing
  the plateau to an experimental binding constant. A rigorous binding free energy needs
  a correction for this; the raw well depth is best used to **compare** systems
  (mutants, temperatures, force fields) calculated the same way.
- a **minimum** at short distance — the bound dimer, with the hydrophobic
  cores in contact. The *position* tells you the contact distance; the *depth*
  tells you how stable it is.
- possibly a **barrier** in between — for example from water that has to be
  squeezed out from between the cores, or from side chains that need to
  rearrange.

Limits to keep in mind:

- **One coordinate only.** The PMF tells you about distance. If something
  else slow is happening (the peptides rotate, flip between parallel and
  anti-parallel arrangements, change secondary structure) and is *not*
  sampled within each window, the PMF is incomplete, even if it looks smooth.
  Poor sampling of these hidden coordinates shows up as slow convergence and
  hysteresis.
- **The force field matters.** The PMF is only as good as the energy
  function behind it. CHARMM36m was chosen because it handles disordered
  peptides well, but results should be compared with experiment when possible.
- **The result depends on temperature and composition.** Window temperature
  is set in `pmf_config.sh` (`TEMPERATURE`) and the same value should be given to
  WHAM (`WHAM_TEMP`); the ion concentration also affects the charged residues.

---

## 8. Putting it together: a checklist

1. Decide the system, the reaction coordinate and the force field.
2. Run **01**: build, pull, extract, choose windows. Look at `selected_configs.dat`:
   do the distances cover the range you want, and are they evenly spaced?
3. Run **02** with a short test (few windows, short times) to catch mistakes,
   then with real settings.
4. Run **04** early, on partial data, to look at the histograms and the shape
   of the PMF.
5. Run **03** to extend, and repeat 04. Stop when the PMF no longer changes
   and the error bars are acceptably small.
6. Report the PMF together with: the reaction coordinate definition, the
   force constant and window spacing, the simulation length per window, the
   time discarded, and the error estimate.

---

## Glossary

| Term | Meaning |
|---|---|
| **Reaction coordinate** | The single quantity (here a distance) along which we compute the free energy |
| **PMF** | Potential of mean force: free energy as a function of the reaction coordinate |
| **Umbrella window** | One simulation with a spring holding the coordinate near `d₀` |
| **Bias potential** | The artificial spring `½ k (d − d₀)²` added to a window |
| **Force constant `k`** | Stiffness of the spring (kJ mol⁻¹ nm⁻²); stiffer = narrower histogram |
| **Steered MD** | Dragging the system along the coordinate to generate starting structures |
| **WHAM** | Weighted Histogram Analysis Method: removes the bias and merges the windows |
| **Overlap** | Neighbouring histograms sharing some distances; required for WHAM to work |
| **Bootstrap** | Repeating the analysis on resampled data to estimate the error |
