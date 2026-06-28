import re
import pandas as pd
import glob
import os

data_folder = os.getenv("ARHG_CSV_DIR", os.path.join(os.path.expanduser("~"), "data", "csvs"))
out_csv     = "all_samples_clean.csv"

# List CSV files inside the folder
csv_files = glob.glob(os.path.join(data_folder, "*.csv"))
print(f"Found {len(csv_files)} CSVs")
for f in csv_files:
  print(" -", os.path.basename(f))

# Columns that are shared (not sample-specific)
BASE_COLS = [
  'Chr', 'Start', 'Stop', 'Ref', 'Alt',
  'Filter', 'Gene', 'MANE HGVS', 'Loc In Gene', 'Effect'
]

# Regex patterns to discover sample IDs from column names
RE_ALT   = re.compile(r"^(?P<sample>.+?)\s+Alt Percentage$")
RE_AD    = re.compile(r"^(?P<sample>.+?)\s+AD Total$")
RE_TNVAF = re.compile(r"^T-N VAF\s*\((?P<sample>.+?)\)$")

def find_samples(columns):
  samples = set()
for col in columns:
  for pat in (RE_ALT, RE_AD, RE_TNVAF):
  m = pat.match(col)
if m:
  samples.add(m.group("sample"))
break
return samples

def build_per_sample(df, sample_id):
  # Expected sample-specific columns
  alt_col = f"{sample_id} Alt Percentage"
ad_col  = f"{sample_id} AD Total"
vaf_col = f"T-N VAF ({sample_id})"

missing = [c for c in (alt_col, ad_col, vaf_col) if c not in df.columns]
if missing:
  # Skip if missing any of the trio for this sample in this file
  return None

# Keep only available base columns (in case some files miss a few)
base_present = [c for c in BASE_COLS if c in df.columns]

# Build subset
sub = df[base_present + [alt_col, ad_col, vaf_col]].copy()
sub = sub.rename(columns={
  alt_col: 'Alt Percentage',
  ad_col:  'AD Total',
  vaf_col: 'T-N VAF'
})
sub['Sample ID'] = sample_id

# Move "Sample ID" to front
ordered_cols = ['Sample ID'] + [c for c in sub.columns if c != 'Sample ID']
return sub[ordered_cols]

# --- Main pipeline over all CSVs ---
all_rows = []
files_processed = 0
samples_made = 0

for filepath in csv_files:
  df = pd.read_csv(filepath)
files_processed += 1

# Filter to PASS
if 'Filter' not in df.columns:
  print(f"[WARN] Skipping {os.path.basename(filepath)}: no 'Filter' column")
continue
df_pass = df.loc[df['Filter'] == 'PASS'].copy()

# Detect sample IDs present in this file
samples = find_samples(df_pass.columns)
if not samples:
  print(f"[INFO] No sample-specific columns detected in {os.path.basename(filepath)}")
continue

for sample_id in samples:
  sub = build_per_sample(df_pass, sample_id)
if sub is not None and not sub.empty:
  sub['File'] = os.path.basename(filepath)  # optional provenance
all_rows.append(sub)
samples_made += 1

# Concatenate and save
if all_rows:
  result = pd.concat(all_rows, ignore_index=True)
# Optional: enforce column order with Sample ID first
cols = ['Sample ID'] + [c for c in result.columns if c != 'Sample ID']
result = result[cols]
result.to_csv(out_csv, index=False)
print(f"Saved {len(result)} rows from {files_processed} CSV(s), {samples_made} sample tables → {out_csv}")
else:
  print(f"No rows found across {files_processed} CSV(s). Check Filter=='PASS' and column name patterns.")
