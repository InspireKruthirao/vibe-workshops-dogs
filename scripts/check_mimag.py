"""Pin down the exact MIMAG rule and the '0' sentinel in the GTDB ref column."""
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1])
hq = df[df.Quality == "high-quality"]
fails = hq[hq.MIMAG == "No"]

print("high-quality rows flagged MIMAG=No:", len(fails))
cols = ["16S total", "23S total", "5S total", "Unique tRNAs", "Total tRNAs",
        "%16S partial", "%23S partial", "%5S partial"]
print(fails[cols].to_string(index=False))

print("\namong the 14 failures, how many have each deficiency:")
print("  16S == 0      :", int((fails["16S total"] == 0).sum()))
print("  23S == 0      :", int((fails["23S total"] == 0).sum()))
print("  5S  == 0      :", int((fails["5S total"] == 0).sum()))
print("  uniq tRNA < 18:", int((fails["Unique tRNAs"] < 18).sum()))

passes = hq[hq.MIMAG == "Yes"]
print("\namong the 1930 MIMAG=Yes rows, any deficiency present?")
print("  16S == 0      :", int((passes["16S total"] == 0).sum()))
print("  23S == 0      :", int((passes["23S total"] == 0).sum()))
print("  5S  == 0      :", int((passes["5S total"] == 0).sum()))
print("  uniq tRNA < 18:", int((passes["Unique tRNAs"] < 18).sum()))

# candidate rule: HQ + 16S & 23S & 5S present + >=18 unique tRNA, but 5S ignored?
for name, rule in {
    "HQ & 16S>0 & 23S>0 & 5S>0 & tRNA>=18":
        (df["16S total"] > 0) & (df["23S total"] > 0) & (df["5S total"] > 0) & (df["Unique tRNAs"] >= 18),
    "HQ & 16S>0 & 23S>0 & tRNA>=18 (5S ignored)":
        (df["16S total"] > 0) & (df["23S total"] > 0) & (df["Unique tRNAs"] >= 18),
    "HQ & tRNA>=18 only":
        df["Unique tRNAs"] >= 18,
    "HQ only":
        pd.Series(True, index=df.index),
}.items():
    pred = (df.Quality == "high-quality") & rule
    print(f"  mismatches for [{name}]: {int(((df.MIMAG == 'Yes') != pred).sum())}")

print("\n--- '0' sentinel in GTDBtk fastani Ref ---")
zero = df["GTDBtk fastani Ref"].astype(str) == "0"
print("rows with ref == '0':", int(zero.sum()))
species = df.Classification.str.split(";").str[6].str[3:]
print("rows with empty species:", int((species == "").sum()))
print("ref=='0' exactly coincides with empty species:",
      bool((zero == (species == "")).all()))
print("\nexample unresolved rows:")
print(df.loc[zero, ["Bin ID", "GTDBtk fastani Ref", "Classification"]].head(3).to_string(index=False))
