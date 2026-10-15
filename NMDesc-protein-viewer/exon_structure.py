"""Map an Ensembl transcript's exons onto protein (residue) coordinates.

The viewer's x-axis is residue position, so an exon track has to be expressed in
residues. Only the coding part of each exon counts: UTR-only exons carry no
residue and are reported separately.

Method, per transcript:
  1. order exons 5'->3' (ascending genomic start on + strand, descending on -)
  2. intersect each exon with the translated span to get its coding length
  3. accumulate coding length to get each exon's CDS interval, 1-based
  4. residues = ceil(cds_start / 3) .. ceil(cds_end / 3)

An exon boundary generally falls inside a codon, so consecutive exons can share
one residue; that shared codon is the real split codon and is kept, not hidden.

The LAST exon-exon junction is reported separately: the 50-nucleotide rule that
defines NMD escape is stated relative to it, so it is what makes a PTC escape or
not. It is given both in CDS nucleotides and in residues.
"""
import math, json, os, urllib.request, urllib.error

ENSEMBL = "https://rest.ensembl.org"


def fetch_transcript(tid, cache_dir="data/cache", retries=3):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"ens_{tid}.json")
    if os.path.exists(path):
        try:
            return json.load(open(path))
        except Exception:
            os.remove(path)
    url = f"{ENSEMBL}/lookup/id/{tid}?expand=1"
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            json.loads(body)
            tmp = path + ".part"
            open(tmp, "wb").write(body)
            os.replace(tmp, path)
            return json.loads(body)
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            if k == retries - 1:
                raise
        except Exception:
            if k == retries - 1:
                raise
    return None


def exon_structure(rec):
    """-> {exons, n_exons, cds_length, last_junction} or None."""
    if not rec or "Exon" not in rec or not rec.get("Translation"):
        return None
    strand = rec.get("strand", 1)
    tr = rec["Translation"]
    cds_lo, cds_hi = min(tr["start"], tr["end"]), max(tr["start"], tr["end"])
    exons = sorted(rec["Exon"], key=lambda e: e["start"], reverse=(strand == -1))

    out, offset = [], 0
    for i, e in enumerate(exons, 1):
        lo, hi = max(e["start"], cds_lo), min(e["end"], cds_hi)
        coding = max(0, hi - lo + 1)
        row = dict(rank=i, exon_id=e["id"], genomic_start=e["start"], genomic_end=e["end"],
                   coding_nt=coding)
        if coding:
            row["cds_start"] = offset + 1
            row["cds_end"] = offset + coding
            row["start"] = int(math.ceil(row["cds_start"] / 3))
            row["end"] = int(math.ceil(row["cds_end"] / 3))
            offset += coding
        out.append(row)

    coding_exons = [e for e in out if e["coding_nt"]]
    if not coding_exons:
        return None
    last_junction = None
    if len(coding_exons) > 1:
        prev = coding_exons[-2]
        last_junction = dict(cds_nt=prev["cds_end"], residue=prev["end"],
                             after_exon_rank=prev["rank"])
    return dict(exons=out, n_exons=len(out), n_coding_exons=len(coding_exons),
                cds_length=offset, last_junction=last_junction, strand=strand)


def check(struct, expected_cds, expected_aa):
    """Agreement between the Ensembl structure and the project / UniProt lengths."""
    if struct is None:
        return "no_structure"
    notes = []
    # the project's cds_length includes the stop codon; Ensembl's translated span does not
    if expected_cds is not None and struct["cds_length"] not in (expected_cds, expected_cds - 3):
        notes.append(f"cds {struct['cds_length']} vs table {expected_cds}")
    aa = struct["cds_length"] // 3
    if expected_aa is not None and abs(aa - expected_aa) > 1:
        notes.append(f"aa {aa} vs uniprot {expected_aa}")
    return "; ".join(notes) if notes else "ok"
