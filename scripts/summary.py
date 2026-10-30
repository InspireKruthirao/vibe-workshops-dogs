"""Final descriptive numbers for data.md."""
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1])
ranks = df.Classification.str.split(";")

print("rows:", len(df), "cols:", df.shape[1])
print("\nSample ids: n =", df.Sample.nunique())
present = set(df.Sample)
missing = [f"D{i:03d}" for i in range(53) if f"D{i:03d}" not in present]
print("gaps in D000-D052:", missing)
vc = df.Sample.value_counts()
print("MAGs per sample: min=%d median=%.0f max=%d" % (vc.min(), vc.median(), vc.max()))

print("\nQuality:", df.Quality.value_counts().to_dict())
print("MIMAG:", df.MIMAG.value_counts().to_dict())
print("Representative:", df.Representative.value_counts().to_dict())
print("MIMAG=Yes among representatives:",
      df[df.Representative == "Yes"].MIMAG.value_counts().to_dict())

print("\nspecies unassigned (s__ empty): %d (%.1f%%)"
      % ((ranks.str[6].str[3:] == "").sum(), (ranks.str[6].str[3:] == "").mean() * 100))
print("genus unassigned: %d" % (ranks.str[5].str[3:] == "").sum())
print("unique species (named):", ranks.str[6].str[3:].replace("", None).nunique())
print("unique genera (named):", ranks.str[5].str[3:].replace("", None).nunique())
print("\nphylum counts:")
print(ranks.str[1].str[3:].value_counts().to_string())

print("\nGenome Size: min=%d median=%d max=%d" %
      (df["Genome Size"].min(), df["Genome Size"].median(), df["Genome Size"].max()))
print("Nr contigs: min=%d median=%d max=%d  (==1 in %d rows)" %
      (df["Nr contigs"].min(), df["Nr contigs"].median(), df["Nr contigs"].max(),
       int((df["Nr contigs"] == 1).sum())))
print("Total tRNAs: min=%d median=%d max=%d" %
      (df["Total tRNAs"].min(), df["Total tRNAs"].median(), df["Total tRNAs"].max()))

print("\nassembly prefix in Original ID:")
print(df["Original ID"].str.split("_SemiBin_").str[0].value_counts().to_string())
print("\nBin ID numeric range:",
      df["Bin ID"].str.extract(r"SHD1_(\d+)")[0].astype(int).agg(["min", "max"]).to_dict())
print("Bin ID unique:", df["Bin ID"].nunique(), "| Original ID unique:", df["Original ID"].nunique())
