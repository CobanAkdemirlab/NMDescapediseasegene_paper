#!/usr/bin/env bash
# Build the input archive the Colab notebook downloads.
# Reads colab_input_inventory.csv (basename,path,bytes,...) and copies each file it
# found on this machine into one flat directory, then tars it.
#
#   bash make_data_bundle.sh [max_mb]
#
# max_mb (default 25) caps per-file size. Files above the cap are listed in
# oversized.txt so they can be deposited separately.

set -euo pipefail
MAX_MB="${1:-25}"
INV="colab_input_inventory.csv"
OUT="nmdesc_colab_data"
[ -f "$INV" ] || { echo "missing $INV"; exit 1; }

rm -rf "$OUT" && mkdir -p "$OUT"
: > oversized.txt

python3 - "$INV" "$OUT" "$MAX_MB" <<'PY'
import csv, os, shutil, sys
inv, out, max_mb = sys.argv[1], sys.argv[2], float(sys.argv[3])
cap = max_mb * 1e6
copied = skipped = big = 0
with open(inv, newline='') as fh, open('oversized.txt', 'w') as ov:
    for r in csv.DictReader(fh):
        p = r.get('path') or ''
        if not p or not os.path.isfile(p):
            skipped += 1; continue
        sz = os.path.getsize(p)
        if sz > cap:
            ov.write(f"{sz}\t{p}\n"); big += 1; continue
        dst = os.path.join(out, os.path.basename(p))
        if not os.path.exists(dst):
            shutil.copy2(p, dst); copied += 1
print(f"copied {copied} | not found on this machine {skipped} | over {max_mb:g} MB {big}")
PY

tar czf "${OUT}.tar.gz" "$OUT"
echo "archive: ${OUT}.tar.gz  $(du -h "${OUT}.tar.gz" | cut -f1)"
echo "files in archive: $(tar tzf "${OUT}.tar.gz" | grep -c -v '/$' || true)"
echo "sha256: $(shasum -a 256 "${OUT}.tar.gz" | cut -d' ' -f1)"
[ -s oversized.txt ] && echo "over cap, deposit separately: $(wc -l < oversized.txt) files (see oversized.txt)"
