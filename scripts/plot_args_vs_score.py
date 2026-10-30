"""Scatter: ARGs carried per MAG against CheckM2 quality score.

    python3 scripts/plot_args_vs_score.py <mimag.csv> <args.csv.gz> <light.png> <dark.png> <table.csv>

The ARG table is RGI/CARD output with one row per detected ARG ORF. It is joined
on `Original ID`, which matches 1,591/1,591 - the `Bin ID` columns differ by a
`.fna.gz` suffix and will NOT join directly. The join must be a LEFT join from
the MAG table: 1,085 MAGs carry no ARG at all and so are simply absent from the
ARG file. Inner-joining silently drops 41% of the MAGs, all of them zeros, and
biases every count upward.

The y axis is symlog with a linear region below 1, because the counts are heavily
tailed: 97% of MAGs carry 5 ARGs or fewer, but the tail runs to 61. A linear axis
would squeeze almost all of the data into the bottom twelfth of the panel, and a
log axis cannot show the 1,085 zeros at all.

Two series (per-MAG dots, binned mean), so slots 1-2 and a legend. Validated with
the dataviz validator under --pairs all in both modes (scatter needs all-pairs):
worst pair CVD dE 24.7 light / 26.8 dark, normal-vision 33.6 / 31.8, all PASS.
"""
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MIMAG, ARGS, OUT_LIGHT, OUT_DARK, OUT_TABLE = sys.argv[1:6]

THEME = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", dots="#2a78d6", mean="#eb6834"),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", dots="#3987e5", mean="#d95926"),
}

DOT_PX = 8.0       # marker diameter - the >= 8px floor for markers
RING_PX = 1.2      # surface ring, so overlapping dots stay countable
LINE_PX = 2.0      # line spec
END_PX = 9.0       # end marker, >= 8px
DPI = 200
PT = 72 / DPI      # one device pixel in points
SCORE_BINS = [0, 50, 60, 70, 80, 90, 95, 100]

mags = pd.read_csv(MIMAG)
args = pd.read_csv(ARGS, usecols=["ORF_ID", "Best_Hit_ARO", "Original ID",
                                  "Cut_Off"], low_memory=False)

per_mag = args.groupby("Original ID").agg(
    args=("ORF_ID", "size"), distinct_aros=("Best_Hit_ARO", "nunique"))

d = (mags.assign(score=mags.Completeness - 5 * mags.Contamination)
         .merge(per_mag, on="Original ID", how="left")
         .fillna({"args": 0, "distinct_aros": 0}))
d["args"] = d["args"].astype(int)
d["distinct_aros"] = d["distinct_aros"].astype(int)
assert len(d) == len(mags), "left join changed the row count"

d[["Bin ID", "Original ID", "Sample", "Completeness", "Contamination", "score",
   "Quality", "args", "distinct_aros", "Classification"]].to_csv(OUT_TABLE,
                                                                 index=False)

binned = d.groupby(pd.cut(d.score, SCORE_BINS, include_lowest=True),
                   observed=True).agg(mags=("args", "size"),
                                      mean=("args", "mean"),
                                      carrying=("args", lambda s: 100 * (s > 0).mean()))
binned["mid"] = [iv.mid for iv in binned.index]

zero = int((d.args == 0).sum())
r_p = np.corrcoef(d.score, d.args)[0, 1]
r_s = np.corrcoef(d.score.rank(), d.args.rank())[0, 1]
worst = d.args.idxmax()


def render(mode, outfile):
    c = THEME[mode]
    fig, ax = plt.subplots(figsize=(10.5, 6.8), dpi=DPI)
    fig.patch.set_facecolor(c["surface"])
    ax.set_facecolor(c["surface"])

    ax.set_yscale("symlog", linthresh=1, linscale=0.32)
    ax.set_ylim(-0.25, 90)
    ax.set_yticks([0, 1, 2, 3, 5, 10, 20, 50])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlim(0, 104)
    ax.set_xticks(range(0, 101, 20))
    ax.tick_params(labelsize=8, colors=c["muted"], length=0, pad=4)
    ax.set_xlabel("Quality score   (completeness - 5 x contamination)",
                  fontsize=9, color=c["ink2"], labelpad=8)
    ax.set_ylabel("ARGs detected in the MAG", fontsize=9, color=c["ink2"],
                  labelpad=8)

    ax.grid(True, color=c["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["axis"])
        ax.spines[side].set_linewidth(0.8)

    ax.scatter(d.score, d.args, s=(DOT_PX * PT) ** 2, facecolor=c["dots"],
               edgecolor=c["surface"], linewidth=RING_PX * PT, zorder=3,
               label=f"One MAG  ({len(d):,})")
    ax.plot(binned.mid, binned["mean"], color=c["mean"], lw=LINE_PX * PT,
            solid_joinstyle="round", solid_capstyle="round", zorder=5,
            marker="o", markersize=END_PX * PT,
            markeredgecolor=c["surface"], markeredgewidth=RING_PX * PT,
            label="Mean per score band")

    leg = ax.legend(loc="upper left", frameon=False, fontsize=8.5,
                    handlelength=1.1, borderpad=0, labelspacing=0.6)
    for txt in leg.get_texts():
        txt.set_color(c["ink2"])

    # Two direct labels only: the zero row, which is the dominant feature and
    # easy to misread as an axis, and the single worst outlier.
    ax.annotate(f"{zero:,} MAGs ({100 * zero / len(d):.0f}%) carry no ARG at all",
                (1.5, 0), xytext=(14, 0.22), fontsize=7.5, color=c["ink2"],
                ha="left", va="center", zorder=6)
    ax.annotate(f"{d.args[worst]} ARGs - {d.Quality[worst]},\n"
                f"score {d.score[worst]:.0f}",
                (d.score[worst], d.args[worst]),
                xytext=(d.score[worst] - 6, d.args[worst] + 18),
                fontsize=7.5, color=c["ink2"], ha="right", va="center",
                linespacing=1.5, zorder=6,
                arrowprops=dict(arrowstyle="-", color=c["axis"], lw=0.7,
                                shrinkA=3, shrinkB=5))

    fig.text(0.052, 0.957, "Quality score predicts ARG detection, not ARG load",
             fontsize=13.5, color=c["ink"], weight="bold", ha="left", va="top")
    fig.text(0.052, 0.912, "\n".join(textwrap.wrap(
             f"ARGs per MAG against quality score, all {len(d):,} dog gut MAGs "
             f"including the {zero:,} that "
             f"carry none. The share carrying at least one ARG climbs steadily "
             f"with score - {binned.carrying.iloc[0]:.0f}% in the lowest band "
             f"to {binned.carrying.iloc[-1]:.0f}% in the highest - but among "
             f"the MAGs that do carry ARGs the median count stays at 1-2 "
             f"throughout. Pearson r = {r_p:.2f}, Spearman {r_s:.2f}.", 150)),
             fontsize=8.5, color=c["ink2"], ha="left", va="top", linespacing=1.6)
    fig.subplots_adjust(left=0.072, right=0.985, top=0.805, bottom=0.112)

    fig.text(0.052, 0.022,
             "Sources: SHD_bins_MIMAG_report.csv + Zenodo 16356977 (RGI/CARD, "
             "all calls Perfect or Strict), left-joined on Original ID  ·  "
             "symlog y below 1 keeps the zeros visible  ·  per-MAG counts in "
             f"{OUT_TABLE.rsplit('/', 1)[-1]}",
             fontsize=7, color=c["muted"], ha="left", va="bottom")

    fig.savefig(outfile, dpi=DPI, facecolor=c["surface"])
    plt.close(fig)


render("light", OUT_LIGHT)
render("dark", OUT_DARK)
print(f"MAGs={len(d)}  with ARGs={int((d.args > 0).sum())}  zero={zero}  "
      f"max={int(d.args.max())}  ARG rows={len(args)}  "
      f"cut-offs={args.Cut_Off.value_counts().to_dict()}")
print(f"pearson={r_p:.3f}  spearman={r_s:.3f}")
print(binned.round(2).to_string())
