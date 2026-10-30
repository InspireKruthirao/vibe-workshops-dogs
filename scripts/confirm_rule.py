"""Confirm the exact rule reproducing the MIMAG column."""
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1])
for r in ("16S", "23S", "5S"):
    df[f"{r} complete"] = df[f"{r} total"] * (100 - df[f"%{r} partial"]) / 100

hq = (df.Completeness > 90) & (df.Contamination < 5)
yes = df.MIMAG == "Yes"

cands = {
    "HQ & 16S total>0 & 23S total>0 & uniq tRNA>=20":
        hq & (df["16S total"] > 0) & (df["23S total"] > 0) & (df["Unique tRNAs"] >= 20),
    "HQ & 16S complete>0 & 23S complete>0 & uniq tRNA>=20":
        hq & (df["16S complete"] > 0) & (df["23S complete"] > 0) & (df["Unique tRNAs"] >= 20),
    "HQ & 16S/23S/5S total>0 & uniq tRNA>=20":
        hq & (df["16S total"] > 0) & (df["23S total"] > 0) & (df["5S total"] > 0) & (df["Unique tRNAs"] >= 20),
    "HQ & uniq tRNA>=20 (no rRNA term)":
        hq & (df["Unique tRNAs"] >= 20),
}
for name, pred in cands.items():
    print(f"mismatches={int((yes != pred).sum()):5d}  {name}")

best = hq & (df["16S total"] > 0) & (df["23S total"] > 0) & (df["Unique tRNAs"] >= 20)
print("\nEXACT MATCH:", bool((yes == best).all()))
print("\nsanity on MIMAG=Yes subset:")
y = df[yes]
print("  min 16S total:", y["16S total"].min(), "| min 23S total:", y["23S total"].min())
print("  min unique tRNAs:", int(y["Unique tRNAs"].min()))
print("  rows with 5S total==0:", int((y["5S total"] == 0).sum()))
print("  rows where 16S all partial (complete==0):", int((y["16S complete"] == 0).sum()))
print("  rows where 23S all partial (complete==0):", int((y["23S complete"] == 0).sum()))
