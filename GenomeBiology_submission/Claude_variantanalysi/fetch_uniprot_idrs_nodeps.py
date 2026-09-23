#!/usr/bin/env python3
"""
Same as fetch_uniprot_idrs.py but uses ONLY the Python standard library
(no pip, no pandas, no requests).

    python3 fetch_uniprot_idrs_nodeps.py variant_all5_0918_1.csv uniprot_idrs.csv [min_len] [extra_accessions]

Writes ID,Start,End (UniProt "Disordered" regions, 1-based inclusive,
>= min_len residues; default 20) AND uniprot_idrs.fasta with each protein's UniProt
canonical sequence, so recompute_idr_columns.py can map UniProt coordinates onto
the Ensembl protein when the two isoforms differ (e.g. RARB 455 vs 448 aa).
extra_accessions: comma-separated accessions missing from the variant list,
e.g. P10070,O95319 (GLI2, CELF2).
"""
import csv, json, sys, time, urllib.request, urllib.error

SRC = sys.argv[1] if len(sys.argv) > 1 else 'variant_all5_0918_1.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'uniprot_idrs.csv'
MIN_LEN = int(sys.argv[3]) if len(sys.argv) > 3 else 20
EXTRA = [a.strip() for a in sys.argv[4].split(',')] if len(sys.argv) > 4 else []

csv.field_size_limit(10**9)          # the variant file has very long CDS columns
with open(SRC, newline='') as f:
    reader = csv.DictReader(f)
    col = 'uniprotswissprot' if 'uniprotswissprot' in reader.fieldnames else 'uniprot'
    accs = sorted({row[col].split('.')[0].strip() for row in reader
                   if row.get(col) and row[col].strip() and row[col].strip().lower() != 'nan'})
accs = sorted(set(accs) | set(EXTRA))
print('%d UniProt accessions' % len(accs))

rows, failed, seqs = [], [], {}
for i, acc in enumerate(accs, 1):
    url = 'https://rest.uniprot.org/uniprotkb/%s.json?fields=ft_region,sequence' % acc
    data = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)
            break
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            time.sleep(2 ** attempt)
    if data is None:
        failed.append(acc)
        continue
    seqs[acc] = data.get('sequence', {}).get('value', '')
    for f in data.get('features', []):
        if f.get('type') != 'Region' or 'disordered' not in f.get('description', '').lower():
            continue
        s = f['location']['start'].get('value')
        e = f['location']['end'].get('value')
        if s is not None and e is not None and e - s + 1 >= MIN_LEN:
            rows.append((acc, int(s), int(e)))
    if i % 25 == 0:
        print('  %d/%d' % (i, len(accs)))
    time.sleep(0.2)

with open(OUT, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['ID', 'Start', 'End'])
    w.writerows(rows)
with open(OUT.rsplit('.', 1)[0] + '.fasta', 'w') as f:
    for a, q in seqs.items():
        f.write('>%s\n%s\n' % (a, q))
print('wrote %s: %d IDRs (>= %d aa) in %d proteins' % (OUT, len(rows), MIN_LEN, len({r[0] for r in rows})))
if failed:
    print('FAILED (re-run to retry):', ', '.join(failed))
