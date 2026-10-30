# Dog gut MAG dataset — structure and first look

Working notes on the dog metagenome-assembled genome (MAG) dataset from the
[vibe-workshops](https://github.com/luispedro/vibe-workshops) intro session.

The dataset arrived with only a short email for documentation, so the first job
was working out the schema. [`data.md`](data.md) is the result: the full column
reference, the rules behind the two derived columns, and the gotchas that bite
if you take the table at face value.

## Contents

| Path | What it is |
|---|---|
| [`data.md`](data.md) | **Start here** — data dictionary for the MAG table |
| `data/SHD_bins_MIMAG_report.csv` | 2,676 MAGs × 20 columns (CheckM2 stats + GTDB-Tk taxonomy) |
| `data/dataset_email.txt` | The original note describing the dataset |
| `figures/` | Rendered charts, light and dark, plus the underlying counts |
| `scripts/` | Profiling and plotting scripts — every number in `data.md` is reproducible |

## Two things worth knowing up front

- **The `MIMAG` column uses non-standard thresholds.** It requires ≥20 distinct
  tRNAs rather than the published standard's 18, and does not require 5S rRNA at
  all. Use the column as given rather than recomputing it from the standard.
- **`GTDBtk fastani Ref` uses the string `0` to mean "no match"**, not an empty
  cell — 470 rows, exactly the rows whose species is unassigned.

Full details and six more gotchas in [`data.md`](data.md).

## Reproducing

Needs `pandas` and `matplotlib`.

```sh
# column profile and the checks behind data.md
python3 scripts/profile_data.py   data/SHD_bins_MIMAG_report.csv
python3 scripts/check_relations.py data/SHD_bins_MIMAG_report.csv
python3 scripts/check_mimag.py    data/SHD_bins_MIMAG_report.csv
python3 scripts/confirm_rule.py   data/SHD_bins_MIMAG_report.csv
python3 scripts/summary.py        data/SHD_bins_MIMAG_report.csv

# the figure
python3 scripts/plot_quality_by_sample.py \
    data/SHD_bins_MIMAG_report.csv \
    figures/quality_by_sample.png \
    figures/quality_by_sample_dark.png \
    figures/quality_by_sample.csv
```

## MAG recovery per dog

![MAG recovery per dog, split by genome quality](figures/quality_by_sample.png)

Recovery is very uneven — 4 to 108 MAGs per sample, median 46 — but the
high-quality share stays fairly stable (50–100%, median 74%) and barely tracks
total yield (r = −0.27). The low-yield dogs are producing *fewer* genomes, not
systematically worse ones.
