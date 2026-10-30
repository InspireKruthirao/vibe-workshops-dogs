"""Inspect the high-quality MIMAG=No rows that simple rRNA/tRNA rules do not explain."""
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1])
for r in ("16S", "23S", "5S"):
    df[f"{r} complete"] = (df[f"{r} total"] * (100 - df[f"%{r} partial"]) / 100).round(2)

explained = ((df["16S complete"] <= 0) | (df["23S complete"] <= 0)
             | (df["Unique tRNAs"] < 18))
odd = df[(df.Quality == "high-quality") & (df.MIMAG == "No") & ~explained]

print("high-quality + MIMAG=No not explained by 16S/23S/tRNA:", len(odd))
show = ["Bin ID", "Completeness", "Contamination", "Nr contigs", "16S complete",
        "23S complete", "5S complete", "Unique tRNAs", "Total tRNAs", "Representative"]
print(odd[show].to_string(index=False))

print("\n5S requirement clearly not enforced:")
yes = df[df.MIMAG == "Yes"]
print("  MIMAG=Yes rows with 5S total == 0 :", int((yes["5S total"] == 0).sum()))
print("  MIMAG=Yes rows with 5S complete==0:", int((yes["5S complete"] == 0).sum()))

print("\nMIMAG=Yes is a strict subset of high-quality:",
      bool((yes.Quality == "high-quality").all()))
print("completeness range among MIMAG=Yes: %.2f - %.2f"
      % (yes.Completeness.min(), yes.Completeness.max()))
print("contamination range among MIMAG=Yes: %.2f - %.2f"
      % (yes.Contamination.min(), yes.Contamination.max()))
print("unique tRNA min among MIMAG=Yes:", int(yes["Unique tRNAs"].min()))

print("\nbest simple approximation (HQ & 16S complete>0 & 23S complete>0 & tRNA>=18):")
pred = (df.Quality == "high-quality") & ~explained
print("  mismatches:", int(((df.MIMAG == "Yes") != pred).sum()), "of", len(df))
