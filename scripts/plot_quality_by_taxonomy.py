"""Stacked bar chart: MAG quality tiers broken down by a GTDB taxonomic rank.

    python3 scripts/plot_quality_by_taxonomy.py <csv> <rank> <light.png> <dark.png> <table.csv>

Two panels on one ordering: absolute counts (left) and quality share (right).
Counts at every rank are heavily lopsided, so the rare taxa are sub-pixel in the
count panel; the share panel is what makes their composition readable.

Below the rank of phylum the tail is long and a share computed on three MAGs is
noise, so each rank has a documented minimum-n cut (see RANKS) and the chart
shows only taxa that clear it. The CSV always holds every taxon at that rank, so
the omitted ones are reported rather than lost.

Quality is ordinal (high > medium), so it takes the same single-hue blue ramp as
plot_quality_by_sample.py, checked with the dataviz validator in ordinal mode on
both surfaces. The light step sits at 2.06:1, below the 3:1 mark floor, so the
chart ships the required relief channel: direct share labels plus the CSV table.
"""
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

CSV, RANK, OUT_LIGHT, OUT_DARK, OUT_TABLE = sys.argv[1:6]

# GTDB lineages are a fixed seven-field string, so the rank is a column index.
RANKS = {
    "phylum": dict(field=1, prefix="p__", min_n=1),
    "class":  dict(field=2, prefix="c__", min_n=5),
    "order":  dict(field=3, prefix="o__", min_n=5),
    "family": dict(field=4, prefix="f__", min_n=10),
    "genus":  dict(field=5, prefix="g__", min_n=20),
}
if RANK not in RANKS:
    sys.exit(f"rank must be one of {', '.join(RANKS)}")
spec = RANKS[RANK]

THEME = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", high="#1c5cab", medium="#86b6ef",
                  on_high="#ffffff", on_medium="#0b0b0b"),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", high="#9ec5f4", medium="#2a78d6",
                 on_high="#0b0b0b", on_medium="#ffffff"),
}

GAP_PX = 2.0      # surface gap between touching segments
RADIUS_PX = 4.0   # rounded data-end, square at the baseline
MAX_BAR_PX = 24.0
FIG_W = 13.0
CHROME_BOT_IN = 0.78     # x label + source line
SUB_WRAP = 178           # subtitle characters per line at 8.5pt over FIG_W
SUB_LINE_IN = 0.16
TICK_STEPS = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500]

df = pd.read_csv(CSV)
taxon = (df["Classification"].str.split(";", expand=True)[spec["field"]]
           .str.replace(spec["prefix"], "", regex=False)
           .replace("", "(unassigned)"))

allcounts = (df.assign(Taxon=taxon)
               .pivot_table(index="Taxon", columns="Quality", values="Bin ID",
                            aggfunc="count")
               .fillna(0).astype(int)
               .rename(columns={"high-quality": "high", "medium-quality": "medium"}))
allcounts["total"] = allcounts.high + allcounts.medium
allcounts["high_pct"] = 100 * allcounts.high / allcounts.total
allcounts = allcounts.sort_values("total", ascending=False)
allcounts.index.name = RANK.capitalize()
allcounts["plotted"] = allcounts.total >= spec["min_n"]
allcounts.round(1).to_csv(OUT_TABLE)

counts = allcounts[allcounts.plotted]
dropped = allcounts[~allcounts.plotted]

# Plotted top-to-bottom in rank order, so the largest taxon sits at the top.
rows = counts.iloc[::-1]
n = len(rows)
xmax = int(counts.total.max())
step = next(s for s in TICK_STEPS if xmax / s <= 7)

PLURAL = {"phylum": "phyla", "class": "classes", "order": "orders",
          "family": "families", "genus": "genera"}


def subtitle():
    head = (f"{len(df):,} dog gut MAGs across {len(allcounts)} GTDB "
            f"{PLURAL[RANK]}, ranked by total recovered.")
    renorm = (f"The right panel renormalises each {RANK} to show how its "
              "quality splits.")
    if len(dropped) == 0:
        return (f"{head} Counts are dominated by {counts.index[0]}, so "
                f"{renorm[0].lower()}{renorm[1:]}")
    shown = int(counts.total.sum())
    return (f"{head} The {len(dropped)} {PLURAL[RANK]} with fewer than "
            f"{spec['min_n']} MAGs are omitted - a share computed on a handful "
            f"of genomes is noise - leaving {shown:,} MAGs "
            f"({100 * shown / len(df):.0f}%) on the chart. {renorm}")


def widest_label_in(labels, fontsize, dpi=200):
    """Width of the longest rendered label, in inches - measured, not guessed.

    Taxon names run to 23 characters at genus level, so a fixed left margin
    silently clips them; the mark spec forbids a clipped label.
    """
    probe = plt.figure(figsize=(1, 1), dpi=dpi)
    renderer = probe.canvas.get_renderer()
    widest = max(probe.text(0, 0, s, fontsize=fontsize)
                 .get_window_extent(renderer=renderer).width for s in labels)
    plt.close(probe)
    return widest / dpi


# Floor keeps the shallow ranks on the same margin they were tuned at.
LEFT_IN = max(1.456, widest_label_in(rows.index, 8.5) + 0.18)

SUBTITLE = textwrap.wrap(subtitle(), SUB_WRAP)
# One subtitle line is the baseline height; each extra line grows the header.
CHROME_TOP_IN = 1.00 + SUB_LINE_IN * (len(SUBTITLE) - 1)
plot_h = n * (0.33 if n <= 12 else 0.17)
FIG_H = plot_h + CHROME_TOP_IN + CHROME_BOT_IN


def rounded_end(ax, xl, xr, yb, yt, r_px, color):
    """Bar with the right (data) end rounded by r_px, square at the baseline."""
    t = ax.transData
    (x0, y0), (x1, y1) = t.transform([(xl, yb), (xr, yt)])
    r = min(r_px, x1 - x0, (y1 - y0) / 2)
    verts = [(x0, y0), (x1 - r, y0), (x1, y0), (x1, y0 + r),
             (x1, y1 - r), (x1, y1), (x1 - r, y1), (x0, y1), (x0, y0)]
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3,
             Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]
    path = t.inverted().transform_path(Path(verts, codes))
    ax.add_patch(PathPatch(path, facecolor=color, edgecolor="none", lw=0, zorder=3))


def draw_stack(ax, i, hi, med, bar_h, gap, c):
    """One stacked row: high segment from the baseline, medium riding on top."""
    yb, yt = i - bar_h / 2, i + bar_h / 2
    if med > 0:
        if hi > 0:
            ax.add_patch(Rectangle((0, yb), hi, bar_h, facecolor=c["high"],
                                   edgecolor="none", lw=0, zorder=3))
        rounded_end(ax, hi + gap, hi + med, yb, yt, RADIUS_PX, c["medium"])
    elif hi > 0:
        rounded_end(ax, 0, hi, yb, yt, RADIUS_PX, c["high"])


def style(ax, c, xlim, ticks, labels, xlabel, show_names):
    ax.set_xlim(0, xlim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8, color=c["muted"])
    ax.set_yticks(range(n))
    if show_names:
        # The panels share the y axis, so only one of them may set labels -
        # a FixedFormatter on the other would silently win and blank both.
        ax.set_yticklabels(rows.index, fontsize=8.5, color=c["ink2"])
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(length=0, pad=5)
    ax.set_xlabel(xlabel, fontsize=8.5, color=c["ink2"], labelpad=7)
    ax.xaxis.grid(True, color=c["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(c["axis"])
    ax.spines["left"].set_linewidth(0.8)


def render(mode, outfile):
    c = THEME[mode]
    fig, (axc, axs) = plt.subplots(
        1, 2, figsize=(FIG_W, FIG_H), dpi=200, sharey=True,
        gridspec_kw=dict(width_ratios=[1.0, 0.82], wspace=0.06))
    fig.patch.set_facecolor(c["surface"])
    for ax in (axc, axs):
        ax.set_facecolor(c["surface"])

    ticks = range(0, xmax + 1, step)
    style(axc, c, xmax * 1.03, ticks, [f"{v:,}" for v in ticks],
          "MAGs recovered", show_names=True)
    style(axs, c, 118, range(0, 101, 25), [f"{v}%" for v in range(0, 101, 25)],
          "Share of MAGs in each quality tier", show_names=False)

    fig.text(0.052, 1 - 0.243 / FIG_H, f"Genome quality by {RANK}",
             fontsize=13.5, color=c["ink"], weight="bold", ha="left", va="top")
    fig.text(0.052, 1 - 0.540 / FIG_H, "\n".join(SUBTITLE), fontsize=8.5,
             color=c["ink2"], ha="left", va="top", linespacing=1.5)
    fig.subplots_adjust(left=LEFT_IN / FIG_W, right=0.995,
                        top=1 - CHROME_TOP_IN / FIG_H,
                        bottom=CHROME_BOT_IN / FIG_H)
    fig.canvas.draw()

    # Layout is fixed now, so pixel-based specs can be converted to data units.
    bb = axc.get_window_extent()
    px_y = bb.height / n                      # pixels per row slot
    bar_h = min(MAX_BAR_PX, 0.62 * px_y) / px_y
    gap_c = GAP_PX / (bb.width / (xmax * 1.03))
    gap_s = GAP_PX / (axs.get_window_extent().width / 118)

    for i, (_, row) in enumerate(rows.iterrows()):
        hi, med, tot = int(row.high), int(row.medium), int(row.total)
        draw_stack(axc, i, hi, med, bar_h, gap_c, c)
        draw_stack(axs, i, row.high_pct, 100 - row.high_pct, bar_h, gap_s, c)

        # Relief channel for the sub-3:1 light step: the share in words, inside
        # the high segment where it fits, otherwise outside the bar end.
        pct, label = row.high_pct, f"{row.high_pct:.0f}%"
        if pct >= 18:
            axs.text(pct - 2.5, i, label, fontsize=7.5, color=c["on_high"],
                     ha="right", va="center", zorder=4)
        else:
            axs.text(102, i, label, fontsize=7.5, color=c["ink2"],
                     ha="left", va="center", zorder=4)
        axs.text(117, i, f"n={tot:,}", fontsize=7.5, color=c["muted"],
                 ha="right", va="center", zorder=4)

    handles = [Rectangle((0, 0), 1, 1, facecolor=c["high"], edgecolor="none"),
               Rectangle((0, 0), 1, 1, facecolor=c["medium"], edgecolor="none")]
    labels = [f"High-quality  ({counts.high.sum():,})",
              f"Medium-quality  ({counts.medium.sum():,})"]
    leg = axc.legend(handles, labels, loc="lower right", frameon=False,
                     fontsize=8.5, handlelength=0.9, handleheight=0.9,
                     borderpad=0, labelspacing=0.5)
    for txt in leg.get_texts():
        txt.set_color(c["ink2"])

    fig.text(0.052, 0.119 / FIG_H,
             "Source: SHD_bins_MIMAG_report.csv  ·  high-quality = >90% complete "
             f"and <5% contamination  ·  counts for every {RANK}, plotted or not, "
             f"in {OUT_TABLE.rsplit('/', 1)[-1]}",
             fontsize=7, color=c["muted"], ha="left", va="bottom")

    fig.savefig(outfile, dpi=fig.dpi, facecolor=c["surface"])
    plt.close(fig)
    return px_y, min(MAX_BAR_PX, 0.62 * px_y)


px_y, bar_px = render("light", OUT_LIGHT)
render("dark", OUT_DARK)
print(f"rank={RANK} plotted={n}/{len(allcounts)} "
      f"MAGs={int(counts.total.sum())}/{len(df)} xmax={xmax} step={step}")
print(f"figure {FIG_W}x{FIG_H:.2f}in  slot {px_y:.1f}px -> bar {bar_px:.1f}px "
      f"(cap {MAX_BAR_PX:.0f})")
print(counts.drop(columns="plotted").round(1).to_string())
if len(dropped):
    print(f"omitted ({len(dropped)} {PLURAL[RANK]}, {int(dropped.total.sum())} MAGs, "
          f"{int(dropped.high.sum())} high-quality)")
