import pandas as pd, sys

v = pd.read_csv("variants_all0901.csv")
print("columns in variants_all0901.csv:", list(v.columns))

cands = [c for c in v.columns
         if "transcript" in c.lower() or c.lower() in ("ensembl_transcript_id","enst","tx","transcript_id")]
tx = None
for c in cands:
  if v[c].astype(str).str.startswith("ENST").mean() > 0.5:
  tx = c
break
if tx is None:
  sys.exit("No Ensembl-transcript column found. Candidates seen: %s" % cands)
print("using transcript column:", tx)

if "key" not in v.columns:
  sys.exit("No 'key' column found (need chrom:pos|REF|ALT).")

v["transcript"] = v[tx].astype(str).str.split(".").str[0]
keep = ["key", "transcript"] + (["source"] if "source" in v.columns else [])
inp = v[keep].dropna(subset=["key", "transcript"])
inp = inp[inp["transcript"].str.startswith("ENST")]
inp.to_csv("pipeline_input_all.csv", index=False)

print("wrote pipeline_input_all.csv:", len(inp), "rows")
if "source" in inp.columns:
  print(inp["source"].value_counts())
