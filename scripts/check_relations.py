"""Check inter-column relationships and derived-column rules in the MAG report."""
import re
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1])

print("--- Quality vs MIMAG completeness/contamination rules ---")
hq = (df.Completeness > 90) & (df.Contamination < 5)
mq = (df.Completeness >= 50) & (df.Contamination < 10)
print("rows where Quality==high-quality matches (comp>90 & cont<5):",
      ((df.Quality == "high-quality") == hq).all())
print("all rows satisfy medium rule (comp>=50 & cont<10):", mq.all())

print("\n--- MIMAG flag vs full MIMAG high-quality rule ---")
rrna = (df["16S total"] >= 1) & (df["23S total"] >= 1) & (df["5S total"] >= 1)
full = hq & rrna & (df["Unique tRNAs"] >= 18)
print("MIMAG=='Yes' equals (high-quality & 3 rRNAs & >=18 tRNA):",
      ((df.MIMAG == "Yes") == full).all())
print("mismatches:", int(((df.MIMAG == "Yes") != full).sum()))
print(pd.crosstab(df.MIMAG, df.Quality))

print("\n--- are 16S/23S/5S total columns distinct? ---")
print("16S == 23S identical:", df["16S total"].equals(df["23S total"]))
print("16S == 5S identical:", df["16S total"].equals(df["5S total"]))
print("rows where 16S != 23S:", int((df["16S total"] != df["23S total"]).sum()))
print("corr:\n", df[["16S total", "23S total", "5S total"]].corr().round(4))

print("\n--- Representative vs Classification ---")
reps = df[df.Representative == "Yes"]
print("n representatives:", len(reps))
print("unique classifications overall:", df.Classification.nunique())
print("unique classifications among reps:", reps.Classification.nunique())
print("reps with duplicated classification:",
      int(reps.Classification.duplicated().sum()))
print("classifications with >1 rep (top 5):")
print(reps.Classification.value_counts().head(5).to_string())
print("unique GTDB ref among reps:", reps["GTDBtk fastani Ref"].nunique())

print("\n--- ID structure ---")
print("Bin ID prefixes:", df["Bin ID"].str.split("_").str[0].value_counts().to_dict())
print("Bin ID all end .fna.gz:", bool(df["Bin ID"].str.endswith(".fna.gz").all()))
print("Sample range:", df.Sample.min(), "->", df.Sample.max(),
      "| n =", df.Sample.nunique())
print("Sample all match D\\d{3}:", bool(df.Sample.str.fullmatch(r"D\d{3}").all()))
print("MAGs per sample: min=%d max=%d mean=%.1f"
      % (df.Sample.value_counts().min(), df.Sample.value_counts().max(),
         df.Sample.value_counts().mean()))

print("\nOriginal ID: does it end with the Sample?",
      bool(df.apply(lambda r: r["Original ID"].endswith(r["Sample"]), axis=1).all()))
binner = df["Original ID"].str.extract(r"_(SemiBin|metabat2?|MaxBin2?|concoct)_",
                                       flags=re.I)[0]
print("binner token counts:", binner.value_counts(dropna=False).to_dict())
pre = df["Original ID"].str.split("_SemiBin_").str[0]
print("assembly prefix counts:", pre.value_counts().to_dict())

print("\n--- GTDB taxonomy depth ---")
ranks = df.Classification.str.split(";")
print("all 7 ranks:", bool((ranks.str.len() == 7).all()))
print("prefix order ok:", bool(ranks.apply(
    lambda p: [x[:3] for x in p] == ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]).all()))
for i, rank in enumerate(["domain", "phylum", "class", "order", "family", "genus", "species"]):
    vals = ranks.str[i].str[3:]
    print(f"  {rank:8s} unique={vals.nunique():4d} empty={int((vals == '').sum()):4d}")
print("domains:", ranks.str[0].value_counts().to_dict())
print("top phyla:", ranks.str[1].str[3:].value_counts().head(5).to_dict())

print("\n--- GTDB ref accession prefixes ---")
print(df["GTDBtk fastani Ref"].str[:4].value_counts().to_dict())
print("unique refs:", df["GTDBtk fastani Ref"].nunique())
