"""Profile the MAG report CSV so we can document its structure."""
import sys

import pandas as pd

path = sys.argv[1]
df = pd.read_csv(path)

print(f"shape: {df.shape}")
print(f"duplicate rows: {df.duplicated().sum()}")
print()

for col in df.columns:
    s = df[col]
    nun = s.nunique(dropna=True)
    print(f"=== {col!r}")
    print(f"    dtype={s.dtype} missing={s.isna().sum()} unique={nun}")
    if pd.api.types.is_numeric_dtype(s):
        print(f"    min={s.min()} max={s.max()} mean={round(float(s.mean()), 3)} median={s.median()}")
        print(f"    all_integral={bool((s.dropna() % 1 == 0).all())}")
    if nun <= 12:
        print(f"    values={s.value_counts(dropna=False).to_dict()}")
    else:
        print(f"    examples={s.dropna().unique()[:3].tolist()}")
    print()
