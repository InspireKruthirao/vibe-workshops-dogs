"""Stacked bar chart: MAG count per dog sample, split by CheckM2 quality tier.

Quality is an ordinal encoding (high > medium), so it takes a single-hue blue
ramp with monotone lightness rather than categorical hues. Both ramps were
checked with the dataviz validator (ordinal mode, light and dark surfaces).
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

CSV, OUT_LIGHT, OUT_DARK, OUT_TABLE = sys.argv[1:5]

THEME = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", high="#1c5cab", medium="#86b6ef"),
    # Dark mode is stepped for the dark surface, not an inverted flip; the ramp
    # anchor flips so the stronger tier stays the more prominent one.
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", high="#9ec5f4", medium="#2a78d6"),
}

GAP_PX = 2.0      # surface gap between touching segments
RADIUS_PX = 4.0   # rounded data-end, square at the baseline
MAX_BAR_PX = 24.0

df = pd.read_csv(CSV)
counts = (df.pivot_table(index="Sample", columns="Quality", values="Bin ID",
                         aggfunc="count")
            .fillna(0).astype(int)
            .rename(columns={"high-quality": "high", "medium-quality": "medium"}))
counts["total"] = counts.high + counts.medium
counts = counts.sort_values("total", ascending=False)
counts.to_csv(OUT_TABLE)

n = len(counts)
ymax = int(counts.total.max())


def rounded_top(ax, xl, xr, yb, yt, r_px, color):
    """Rectangle with the top corners rounded by r_px, square at the bottom."""
    t = ax.transData
    (x0, y0), (x1, y1) = t.transform([(xl, yb), (xr, yt)])
    r = min(r_px, (x1 - x0) / 2, y1 - y0)
    verts = [(x0, y0), (x0, y1 - r), (x0, y1), (x0 + r, y1),
             (x1 - r, y1), (x1, y1), (x1, y1 - r), (x1, y0), (x0, y0)]
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3,
             Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]
    path = t.inverted().transform_path(Path(verts, codes))
    ax.add_patch(PathPatch(path, facecolor=color, edgecolor="none", lw=0, zorder=3))


def render(mode, outfile):
    c = THEME[mode]
    fig, ax = plt.subplots(figsize=(13, 5.6), dpi=200)
    fig.patch.set_facecolor(c["surface"])
    ax.set_facecolor(c["surface"])

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, ymax * 1.06)
    ax.set_yticks(range(0, ymax + 1, 25))
    ax.set_xticks(range(n))
    ax.set_xticklabels(counts.index, rotation=90, fontsize=6.5, color=c["muted"])
    ax.tick_params(axis="y", labelsize=8, colors=c["muted"], length=0)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_ylabel("MAGs recovered", fontsize=9, color=c["ink2"], labelpad=8)

    ax.yaxis.grid(True, color=c["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(c["axis"])
    ax.spines["bottom"].set_linewidth(0.8)

    fig.text(0.055, 0.955, "MAG recovery per dog, split by genome quality",
             fontsize=13.5, color=c["ink"], weight="bold", ha="left", va="top")
    fig.text(0.055, 0.905,
             f"{len(df):,} metagenome-assembled genomes from {n} faecal samples, "
             "sorted by total recovered. Tiers are the MIMAG completeness/"
             "contamination bands as scored by CheckM2.",
             fontsize=8.5, color=c["ink2"], ha="left", va="top")
    fig.subplots_adjust(left=0.055, right=0.995, top=0.835, bottom=0.145)
    fig.canvas.draw()

    # Layout is fixed now, so pixel-based specs can be converted to data units.
    bb = ax.get_window_extent()
    px_x = bb.width / (n)              # pixels per x slot
    px_y = bb.height / (ymax * 1.06)   # pixels per y unit
    bar_w = min(MAX_BAR_PX, 0.62 * px_x) / px_x
    gap = GAP_PX / px_y

    for i, (_, row) in enumerate(counts.iterrows()):
        xl, xr = i - bar_w / 2, i + bar_w / 2
        hi, med = int(row.high), int(row.medium)
        if med > 0:
            if hi > 0:
                ax.add_patch(Rectangle((xl, 0), bar_w, hi, facecolor=c["high"],
                                       edgecolor="none", lw=0, zorder=3))
            rounded_top(ax, xl, xr, hi + gap, hi + med, RADIUS_PX, c["medium"])
        elif hi > 0:
            rounded_top(ax, xl, xr, 0, hi, RADIUS_PX, c["high"])

    handles = [Rectangle((0, 0), 1, 1, facecolor=c["high"], edgecolor="none"),
               Rectangle((0, 0), 1, 1, facecolor=c["medium"], edgecolor="none")]
    labels = [f"High-quality  ({counts.high.sum():,})",
              f"Medium-quality  ({counts.medium.sum():,})"]
    leg = ax.legend(handles, labels, loc="upper right", frameon=False,
                    fontsize=8.5, handlelength=0.9, handleheight=0.9,
                    borderpad=0, labelspacing=0.5)
    for txt in leg.get_texts():
        txt.set_color(c["ink2"])

    fig.text(0.055, 0.022,
             "Source: SHD_bins_MIMAG_report.csv  ·  high-quality = >90% complete "
             "and <5% contamination  ·  per-sample counts in quality_by_sample.csv",
             fontsize=7, color=c["muted"], ha="left", va="bottom")

    fig.savefig(outfile, dpi=fig.dpi, facecolor=c["surface"])
    plt.close(fig)
    return px_x, min(MAX_BAR_PX, 0.62 * px_x)


px_x, bar_px = render("light", OUT_LIGHT)
render("dark", OUT_DARK)
print(f"slot width {px_x:.1f}px -> bar {bar_px:.1f}px (cap {MAX_BAR_PX:.0f})")
print(f"samples={n} ymax={ymax} total={int(counts.total.sum())}")
print(counts.head(3).to_string())
