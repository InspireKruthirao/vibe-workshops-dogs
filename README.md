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
| `figures/` | Rendered charts, light and dark, plus the underlying counts (`*.csv`) |
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

# the figures
python3 scripts/plot_quality_by_sample.py \
    data/SHD_bins_MIMAG_report.csv \
    figures/quality_by_sample.png \
    figures/quality_by_sample_dark.png \
    figures/quality_by_sample.csv

# one per taxonomic rank: phylum | class | order | family | genus
for rank in phylum family genus; do
    python3 scripts/plot_quality_by_taxonomy.py \
        data/SHD_bins_MIMAG_report.csv "$rank" \
        "figures/quality_by_$rank.png" \
        "figures/quality_by_${rank}_dark.png" \
        "figures/quality_by_$rank.csv"
done
```

## MAG recovery per dog

![MAG recovery per dog, split by genome quality](figures/quality_by_sample.png)

Recovery is very uneven — 4 to 108 MAGs per sample, median 46 — but the
high-quality share stays fairly stable (50–100%, median 74%) and barely tracks
total yield (r = −0.27). The low-yield dogs are producing *fewer* genomes, not
systematically worse ones.

## Genome quality by taxonomy

`scripts/plot_quality_by_taxonomy.py` takes the rank as an argument and works at
any of phylum / class / order / family / genus. Each chart pairs absolute counts
with the renormalised share, because counts at every rank are lopsided enough to
render the rare taxa sub-pixel. Below phylum the tail is long and a share
computed on three genomes is noise, so each rank has a minimum-n cut (family
≥ 10, genus ≥ 20); the CSV always lists **every** taxon with a `plotted` flag,
so the omitted ones are reported rather than lost.

### Phylum

![Genome quality by phylum, counts and share](figures/quality_by_phylum.png)

All 2,676 MAGs are Bacteria over 11 phyla, every row assigned. The counts are
lopsided — Bacillota_A alone is 1,460 of them (55%) — but quality looks flat:
every phylum with ≥ 20 MAGs sits between 62% and 82% high-quality.
Fusobacteriota (62%) and Campylobacterota (66%) are the weak ends, and the
dominant Bacillota_A (70%) is itself below the dataset average of 73%, so the
phyla that bin most abundantly are not the ones that bin best. The three
single-digit phyla should be read as anecdotes.

### Family

![Genome quality by family, counts and share](figures/quality_by_family.png)

23 of 41 families clear the ≥ 10 cut, covering 2,613 MAGs (98%). Lachnospiraceae
is 1,008 MAGs on its own — 38% of the dataset — at 69% high-quality.
**Clostridiaceae (55%, n=56) and Helicobacteraceae (59%, n=32) are the problem
families**; Acidaminococcaceae (96%, n=28), Burkholderiaceae_A and
Erysipelotrichaceae (both 85%) sit at the top.

### Genus

![Genome quality by genus, counts and share](figures/quality_by_genus.png)

40 of 144 genera clear the ≥ 20 cut, covering 2,107 MAGs (79%). This is where the
variation actually lives: **Clostridium is 47% high-quality (n=30) and Roseburia
50% (n=20), against Eisenbergiella at 97% (n=36) and Amedibacillus at 95%
(n=42)**. Note that *Blautia*, the single most abundant genus (n=171), is one of
the worst at 56%, while its sibling Blautia_A manages 62% — the dominant
Lachnospiraceae lineages are the hardest ones to bin cleanly.

### Quality variation is a genus-level property

Holding the comparison to taxa with ≥ 20 MAGs at every rank, the spread in
high-quality share widens sharply as you go down the tree:

| Rank | Taxa (n ≥ 20) | High-quality share | Spread |
|---|---|---|---|
| Phylum | 8 | 62–82% | 19 pts |
| Family | 20 | 55–96% | 41 pts |
| Genus | 40 | 47–100% | 53 pts |

Phylum-level aggregation hides most of it. If you need a quality expectation for
a lineage, take it from the genus table, not the phylum one.
