#!/usr/bin/env python3
"""Plot the WHAM PMF and the umbrella histograms (overlap check).
usage: plot_pmf.py <analysis_dir> [unit]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1]
unit = {"kJ": "kJ/mol", "kCal": "kcal/mol", "kT": "kT"}.get(sys.argv[2] if len(sys.argv) > 2 else "kJ", "kJ/mol")

def read(path):
    rows = [l.split() for l in open(path) if l.strip() and l[0] not in "#@"]
    return np.array(rows, dtype=float)

p = read(os.path.join(d, "profile.xvg"))
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(p[:, 0], p[:, 1], lw=2)
if p.shape[1] > 2:
    ax.fill_between(p[:, 0], p[:, 1] - p[:, 2], p[:, 1] + p[:, 2], alpha=0.3)
ax.set_xlabel("Distance (nm)"); ax.set_ylabel(f"PMF ({unit})")
fig.tight_layout(); fig.savefig(os.path.join(d, "pmf.png"), dpi=200)

h = read(os.path.join(d, "histo.xvg"))
fig, ax = plt.subplots(figsize=(6, 4))
for k in range(1, h.shape[1]):
    ax.plot(h[:, 0], h[:, k], lw=1)
ax.set_xlabel("Distance (nm)"); ax.set_ylabel("Counts")
ax.set_title("Umbrella histograms (neighbours must overlap)")
fig.tight_layout(); fig.savefig(os.path.join(d, "histograms.png"), dpi=200)
