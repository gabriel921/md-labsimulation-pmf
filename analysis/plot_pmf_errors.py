#!/usr/bin/env python3
"""PMF with bootstrap errors.   usage: plot_pmf_errors.py <analysis_dir> [kJ|kCal|kT]
Reads <analysis_dir>/bsResult.xvg (gmx wham -bsres: column 1 = coordinate, 2 = average profile,
3 = standard deviation) and writes <analysis_dir>/pmf_errors.png."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1]
unit = sys.argv[2] if len(sys.argv) > 2 else "kJ"
rows = []
for line in open(os.path.join(d, "bsResult.xvg")):
    if line.startswith(("#", "@")) or not line.strip():
        continue
    rows.append([float(v) for v in line.split()])
a = np.array(rows)
if a.ndim != 2 or a.shape[1] < 3:
    sys.exit("bsResult.xvg has no error column")
x, y, err = a[:, 0], a[:, 1], a[:, 2]

ylabel = {"kJ": "PMF (kJ/mol)", "kCal": "PMF (kcal/mol)", "kT": "PMF (kT)"}.get(unit, "PMF")
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x, y, lw=1.8, label="PMF (bootstrap average)")
ax.fill_between(x, y - err, y + err, alpha=0.3, label="± 1 standard deviation")
ax.set_xlabel("Reaction coordinate (nm)")
ax.set_ylabel(ylabel)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(d, "pmf_errors.png"), dpi=200)
print("max error: %.3f   mean error: %.3f" % (err.max(), err.mean()))
