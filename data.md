# Dog MAG dataset — data structure

Documentation of the files in `data/`, reconstructed by profiling the data itself.
Luis's email gives the provenance but not the schema, so everything below was
verified against all 2,676 rows. Claims that are inferred rather than verifiable
from the files are marked **[inferred]**.

## Source

Downloaded from the workshop repo
[`luispedro/vibe-workshops`](https://github.com/luispedro/vibe-workshops),
path `1__intro/data/`.

From `dataset_email.txt`:

> This is our dog MAG dataset. It has some basic stats + taxonomy for each MAG.
> Stats are from CheckM2 and the classification is from GTDB-Tk.

So: **MAGs** (metagenome-assembled genomes) recovered from dog gut metagenomes.
Completeness/contamination come from **CheckM2**, taxonomy from **GTDB-Tk**.
The `SHD` prefix and the dog subject matter point to the ShanghaiDogs project
(`luispedro/ShanghaiDogs`) **[inferred]**.

## Files

| File | Size | Contents |
|---|---|---|
| `data/SHD_bins_MIMAG_report.csv` | 658 KB | 2,676 rows × 20 columns, plain CSV with a header row |
| `data/dataset_email.txt` | 278 B | Luis's note describing provenance |

No missing values anywhere in the table, no duplicate rows, no quoted fields,
no embedded commas. It loads cleanly with defaults:

```python
import pandas as pd
df = pd.read_csv("data/SHD_bins_MIMAG_report.csv")   # (2676, 20)
```

## Grain and keys

**One row per MAG.** 2,676 MAGs drawn from 52 samples.

Two unique identifier columns, both 2,676 distinct values:

- `Bin ID` — the curated catalogue ID, e.g. `SHD1_2081.fna.gz`. Numbering runs
  `SHD1_0000` … `SHD1_2675`. Includes the `.fna.gz` suffix, so strip it to get
  a bare ID.
- `Original ID` — the pre-curation ID from the assembly/binning pipeline,
  e.g. `10G_20G_SemiBin_21_D000`. Pattern:
  `<assembly>_SemiBin_<bin number>_<sample>`. Every row ends with its own
  `Sample` value, and every row was binned with **SemiBin**.

Assembly prefixes in `Original ID`, which look like read-type strategies **[inferred]**:

| Prefix | Rows | Likely meaning |
|---|---|---|
| `ALL_LR` | 1,569 | long-read assembly |
| `ALL_multi` | 642 | hybrid / multi-read-type |
| `ALL_SR` | 355 | short-read assembly |
| `10G_20G` | 110 | subsampled depth (10 Gbp / 20 Gbp) |

## Column reference

| # | Column | Type | Range / values | Meaning |
|---|---|---|---|---|
| 1 | `Bin ID` | str | 2,676 unique | Catalogue MAG ID, `SHD1_####.fna.gz` |
| 2 | `Representative` | str | `Yes` 320 / `No` 2,356 | MAG is the representative of its species-level cluster **[inferred]** |
| 3 | `Completeness` | float | 50.05 – 100.0 (median 98.40) | CheckM2 estimated completeness, % |
| 4 | `Contamination` | float | 0.0 – 9.83 (median 0.69) | CheckM2 estimated contamination, % |
| 5 | `Genome Size` | int | 582,373 – 6,077,103 (median 2.66 Mbp) | Assembly length in **base pairs** |
| 6 | `Quality` | str | `high-quality` 1,944 / `medium-quality` 732 | Completeness/contamination tier — derived, see below |
| 7 | `Nr contigs` | int | 1 – 382 (median 5) | Contig count; 691 MAGs are single-contig (complete circular genomes) |
| 8 | `Classification` | str | 299 unique | Full 7-rank GTDB-Tk lineage, `;`-separated |
| 9 | `GTDBtk fastani Ref` | str | 225 unique | Closest GTDB reference genome accession — **uses `0` as "no match"**, see gotchas |
| 10 | `16S total` | float | 0 – 17 (median 5) | Number of 16S rRNA genes detected |
| 11 | `23S total` | float | 0 – 17 (median 5) | Number of 23S rRNA genes detected |
| 12 | `5S total` | float | 0 – 16 (median 5) | Number of 5S rRNA genes detected |
| 13 | `%16S partial` | float | 0 – 100 (median 0) | Percentage of the detected 16S copies that are partial/fragmentary |
| 14 | `%23S partial` | float | 0 – 100 (median 0) | Same, for 23S |
| 15 | `%5S partial` | float | 0 – 100 (median 0) | Same, for 5S |
| 16 | `Unique tRNAs` | int | 5 – 25 (median 22) | Distinct tRNA isotypes found (max 25) |
| 17 | `Total tRNAs` | int | 7 – 149 (median 60) | Total tRNA gene copies |
| 18 | `MIMAG` | str | `Yes` 1,930 / `No` 746 | Meets this project's MIMAG high-quality bar — derived, see below |
| 19 | `Sample` | str | `D000` – `D052`, 52 values | Source sample (dog) |
| 20 | `Original ID` | str | 2,676 unique | Pipeline ID, see above |

The three `*S total` columns are stored as floats only because of how they were
written out — every value is a whole number, so they can be cast to int.

## Derived columns

Two columns are computable from the others. Both rules were checked against all
2,676 rows and match exactly.

### `Quality` — standard MIMAG completeness/contamination tiers

```
high-quality    ⟺  Completeness > 90  and  Contamination < 5
medium-quality  ⟺  otherwise (all rows satisfy Completeness >= 50 and Contamination < 10)
```

The table has already been filtered to ≥50% completeness and <10% contamination,
so there are no low-quality MAGs — the `low-quality` tier simply never appears.

### `MIMAG` — this project's high-quality flag

```
MIMAG == "Yes"  ⟺  Completeness > 90
                    and Contamination < 5
                    and 16S total > 0
                    and 23S total > 0
                    and Unique tRNAs >= 20
```

**Exact match, 0 mismatches in 2,676 rows.** Note two deliberate deviations from
the published MIMAG standard, which are easy to trip over:

- **5S rRNA is not required.** The formal standard requires all three rRNA genes;
  here 459 MAGs flagged `MIMAG = Yes` have `5S total == 0`, and 676 have no
  *complete* 5S copy.
- **The tRNA threshold is 20, not 18.** The published standard asks for ≥18
  distinct tRNAs; this column uses ≥20. Six MAGs with 19 unique tRNAs and
  otherwise-passing stats are flagged `No`, which is what gives the rule away.

`MIMAG = Yes` is a strict subset of `Quality = high-quality`: no medium-quality
MAG is ever flagged `Yes`, and 14 high-quality MAGs are flagged `No`.

Because the thresholds are non-standard, **use the `MIMAG` column as given**
rather than recomputing it from a textbook definition.

## Taxonomy format

`Classification` is a GTDB-Tk lineage string with exactly 7 ranks, always in
order and always present:

```
d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola sp900546645
```

Split on `;` and strip the 3-character `d__`/`p__`/… prefixes:

```python
ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
df[ranks] = df.Classification.str.split(";", expand=True).apply(lambda c: c.str[3:])
```

Diversity at each rank:

| Rank | Distinct named | Unassigned rows |
|---|---|---|
| domain | 1 (all `Bacteria`) | 0 |
| phylum | 11 | 0 |
| class | 12 | 0 |
| order | 24 | 0 |
| family | 41 | 0 |
| genus | 143 | 4 |
| species | 224 | **470 (17.6%)** |

Everything is bacterial — no archaea, and no eukaryotes or viruses. Phylum
composition is dominated by Firmicutes-lineage clades under GTDB's split naming:
`Bacillota_A` 1,460, `Bacteroidota` 330, `Bacillota` 309, `Fusobacteriota` 173,
`Pseudomonadota` 156, `Actinomycetota` 114, `Bacillota_C` 82,
`Campylobacterota` 44, then `Bacillota_B` 5, `Desulfobacterota` 2,
`Deferribacterota` 1.

## Gotchas

These are the things the email doesn't mention and that will bite you:

1. **`GTDBtk fastani Ref` uses the string `0` as a sentinel**, not an empty cell.
   470 rows have `0`, meaning GTDB-Tk found no reference within the ANI
   threshold. Those 470 rows are *exactly* the rows with an empty species
   (`s__`) — verified one-to-one. Convert it before use:

   ```python
   df["GTDBtk fastani Ref"] = df["GTDBtk fastani Ref"].replace("0", pd.NA)
   ```

   Also note that pandas may read this column as a string while the `0` makes
   naive type inference tempting — do not let it become numeric.

2. **Unassigned species are empty strings, not NaN.** `s__` with nothing after it
   means "novel / unnamed species". `pd.read_csv` gives you `"s__"`, which
   survives `.isna()` checks silently. 17.6% of MAGs are in this state, so any
   species-level grouping needs them handled explicitly or they collapse into one
   bogus category.

3. **`Classification` is not a unique key for species clusters.** There are 320
   representatives but only 283 distinct lineages among them, because unnamed
   species share the identical `...;s__` string while being genuinely different
   clusters. One lineage
   (`p__Bacillota_A;…;g__CAG-269;s__`) has 10 separate representatives. Cluster
   on `Bin ID` of the representative, not on the lineage string.

4. **`Genome Size` is in base pairs**, not Mbp or kbp. Divide by 1e6 for plots.

5. **Sample IDs have a gap.** 52 samples labelled `D000`–`D052`, but **`D009` is
   absent**. Don't assume the numbering is contiguous or that there are 53 dogs.

6. **MAG counts per sample are very uneven** — min 4, median 46, max 108. Any
   per-sample aggregate needs normalising, and cross-sample comparisons of raw
   counts will be misleading.

7. **The table is redundant across samples.** The same organism recovered from
   different dogs appears as multiple rows. For any analysis that should treat
   each organism once, filter to `Representative == "Yes"` (320 rows) — of which
   258 are `MIMAG == "Yes"`.

8. **`Unique tRNAs` caps at 25**, so it is not a completeness proxy at the top of
   its range; many good MAGs sit at the ceiling.

## Common subsets

```python
reps  = df[df.Representative == "Yes"]      #  320 — one per species cluster
hq    = df[df.Quality == "high-quality"]    # 1944 — comp >90, cont <5
mimag = df[df.MIMAG == "Yes"]               # 1930 — project HQ bar
closed = df[df["Nr contigs"] == 1]          #  691 — single-contig genomes
novel = df[df["GTDBtk fastani Ref"] == "0"] #  470 — no GTDB species match
```

## Companion table: `data/ARGs/SHD1_ARGs_MAGs.csv.gz`

RGI/CARD antimicrobial-resistance gene calls for the same MAGs, from
[Zenodo record 16356977](https://zenodo.org/records/16356977). 2.9 MB gzipped,
13 MB open; 4,702 rows × 34 columns. Original filename inside the archive is
`RGI_ARGs-MAGs.csv`.

**Grain: one row per detected ARG ORF**, not per MAG. 4,702 hits over 1,591
MAGs, 210 distinct `Best_Hit_ARO` values.

Joining it to the MIMAG report has two traps:

1. **Join on `Original ID`, not `Bin ID`.** Both files have a `Bin ID` column,
   but the ARG file drops the suffix: `SHD1_2081` there against
   `SHD1_2081.fna.gz` here. A direct `Bin ID` merge matches **zero** rows and
   fails silently. `Original ID` matches 1,591/1,591 exactly.
2. **It must be a LEFT join from the MAG table.** 1,085 of the 2,676 MAGs (41%)
   carry no ARG and are therefore *absent from the file entirely*. An inner join
   drops them — and because every dropped row is a zero, every count, mean and
   correlation computed afterwards is biased upward.

```python
per_mag = args.groupby("Original ID").agg(args=("ORF_ID", "size"))
d = mags.merge(per_mag, on="Original ID", how="left").fillna({"args": 0})
```

Other notes on the file:

- `Cut_Off` is `Strict` (3,728) or `Perfect` (974). **No `Loose` calls**, so the
  usual RGI confidence filter is already applied upstream — do not re-filter and
  report it as your own QC. 78 rows carry `Nudged == True` (Loose promoted to
  Strict by RGI's own nudging).
- `ORF_ID` is unique across all 4,702 rows and embeds an absolute path from the
  original analysis (`/data/anna/animal_metagenome/long-mg-dog/...`), which is
  the only place the contributor's directory layout appears.
- `Quality` is duplicated from the MIMAG report and agrees with it, so it is safe
  to ignore in favour of the authoritative column.
- Counting rows per MAG and counting *distinct* `Best_Hit_ARO` per MAG differ for
  289 MAGs — those carry the same gene more than once. Pick deliberately: hits
  measure copies, distinct AROs measure repertoire.
- Three heavy columns (`Predicted_DNA`, `Predicted_Protein`,
  `CARD_Protein_Sequence`) hold full sequences and are most of the file size.
  Pass `usecols=` or the read is needlessly slow.

## Open questions for Luis

Things that genuinely cannot be determined from the files:

- Confirm the assembly prefixes in `Original ID` (`ALL_LR`, `ALL_SR`,
  `ALL_multi`, `10G_20G`) — are these read types and subsampled depths?
- Is `Representative` dereplication at the usual 95% ANI species level, and
  which tool produced it?
- Why is `D009` missing — excluded sample, or renumbering?
- Is the `Unique tRNAs >= 20` threshold (rather than the standard 18) intentional?
- Which CheckM2 and GTDB-Tk versions / GTDB release were used? This matters for
  the taxonomy, since GTDB renames clades between releases.
