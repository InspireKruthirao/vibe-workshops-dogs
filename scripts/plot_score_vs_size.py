"""Scatter: CheckM2 quality score against assembled genome size, one dot per MAG.

    python3 scripts/plot_score_vs_size.py <csv> <light.png> <dark.png> <table.csv>

Quality score is the standard CheckM2 penalty, completeness - 5 x contamination.
Both axes are magnitude on one series, so the mark takes categorical slot 1 and
there is no legend box - the title names what is plotted. Slot 1 clears the 3:1
mark-contrast floor on both surfaces (4.30:1 light, 4.79:1 dark), so no relief
channel is obligated; the CSV ships anyway to match the other figures.

Colouring the dots by the Quality column was considered and rejected: that column
is derived from completeness and contamination, so it would spend the identity
channel re-drawing a threshold that the y position already carries.
"""
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSV, OUT_LIGHT, OUT_DARK, OUT_TABLE = sys.argv[1:5]

THEME = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", series="#2a78d6"),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", series="#3987e5"),
}

DOT_PX = 8.0      # marker diameter - the >= 8px floor for markers
RING_PX = 1.2     # surface ring, so overlapping dots stay countable
DPI = 200
PT = 72 / DPI     # one device pixel in points

df = pd.read_csv(CSV)
size = df["Genome Size"] / 1e6
score = df["Completeness"] - 5 * df["Contamination"]

out = pd.DataFrame({"Bin ID": df["Bin ID"], "Sample": df["Sample"],
                    "Genome Size (Mb)": size.round(3),
                    "Completeness": df["Completeness"],
                    "Contamination": df["Contamination"],
                    "Quality score": score.round(2),
                    "Quality": df["Quality"],
                    "Classification": df["Classification"]})
out.to_csv(OUT_TABLE, index=False)

hq = df.Quality == "high-quality"
r_all = np.corrcoef(size, score)[0, 1]
r_hq = np.corrcoef(size[hq], score[hq])[0, 1]
r_mq = np.corrcoef(size[~hq], score[~hq])[0, 1]

# The two extremes get the only direct labels - sparing is what makes them work.
big = size.idxmax()
worst = score.idxmin()


def phylum_of(i):
    return df.Classification[i].split(";")[1].replace("p__", "")


def render(mode, outfile):
    c = THEME[mode]
    fig, ax = plt.subplots(figsize=(10.5, 6.6), dpi=DPI)
    fig.patch.set_facecolor(c["surface"])
    ax.set_facecolor(c["surface"])

    ax.set_xlim(0.3, 6.4)
    ax.set_ylim(0, 103)
    ax.set_xticks(np.arange(1, 7))
    ax.set_xticklabels([f"{v} Mb" for v in range(1, 7)])
    ax.set_yticks(range(0, 101, 20))
    ax.tick_params(labelsize=8, colors=c["muted"], length=0, pad=4)
    ax.set_xlabel("Assembled genome size", fontsize=9, color=c["ink2"], labelpad=8)
    ax.set_ylabel("Quality score   (completeness - 5 x contamination)",
                  fontsize=9, color=c["ink2"], labelpad=8)

    ax.grid(True, color=c["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["axis"])
        ax.spines[side].set_linewidth(0.8)

    ax.scatter(size, score, s=(DOT_PX * PT) ** 2, facecolor=c["series"],
               edgecolor=c["surface"], linewidth=RING_PX * PT, zorder=3)

    # The largest MAG sits inside the saturated top band, so its label goes to
    # clear space with a leader line rather than on top of its neighbours.
    ax.annotate(f"Largest MAG: {size[big]:.2f} Mb\n{phylum_of(big)}, "
                f"score {score[big]:.0f}",
                (size[big], score[big]), xytext=(5.9, 38),
                fontsize=7.5, color=c["ink2"], ha="right", va="center",
                linespacing=1.5, zorder=4,
                arrowprops=dict(arrowstyle="-", color=c["axis"], lw=0.7,
                                shrinkA=2, shrinkB=5))
    ax.annotate(f"Lowest score: {score[worst]:.1f}\n{phylum_of(worst)}, "
                f"{df.Contamination[worst]:.1f}% contamination",
                (size[worst], score[worst]),
                xytext=(size[worst] + 0.12, score[worst]),
                fontsize=7.5, color=c["ink2"], ha="left", va="center",
                linespacing=1.5, zorder=4)

    fig.text(0.052, 0.957, "Genome size barely predicts MAG quality",
             fontsize=13.5, color=c["ink"], weight="bold", ha="left", va="top")
    fig.text(0.052, 0.908, "\n".join(textwrap.wrap(
             f"CheckM2 quality score against assembled size for all {len(df):,} "
             f"dog gut MAGs. The two track weakly across the whole set "
             f"(r = {r_all:.2f}) for a mechanical reason - an incomplete MAG "
             f"assembles shorter, and completeness is most of the score. Once "
             f"completeness is already above 90%, the relationship disappears: "
             f"within the high-quality MAGs r = {r_hq:.2f}.", 150)),
             fontsize=8.5, color=c["ink2"], ha="left", va="top", linespacing=1.6)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.825, bottom=0.115)

    fig.text(0.052, 0.022,
             "Source: SHD_bins_MIMAG_report.csv  ·  score = completeness - 5 x "
             "contamination, so 100 means complete and uncontaminated  ·  "
             "one dot per MAG, values in "
             f"{OUT_TABLE.rsplit('/', 1)[-1]}",
             fontsize=7, color=c["muted"], ha="left", va="bottom")

    fig.savefig(outfile, dpi=DPI, facecolor=c["surface"])
    plt.close(fig)


render("light", OUT_LIGHT)
render("dark", OUT_DARK)
print(f"n={len(df)}  size {size.min():.2f}-{size.max():.2f} Mb  "
      f"score {score.min():.1f}-{score.max():.1f}")
print(f"pearson all={r_all:.3f}  high-quality={r_hq:.3f}  medium-quality={r_mq:.3f}")
print(f"median size {size.median():.2f} Mb (high-quality {size[hq].median():.2f}, "
      f"medium-quality {size[~hq].median():.2f})")
print(f"median score {score.median():.1f} (high-quality {score[hq].median():.1f}, "
      f"medium-quality {score[~hq].median():.1f})")
