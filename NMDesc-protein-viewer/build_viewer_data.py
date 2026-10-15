"""Build the per-gene JSON payloads behind the NMDesc protein viewer.

One record per (gene, design group) of the gene set 4 case list. Every payload is
self-contained so the static page can fetch exactly one file per selection.

Sources, all public and queried once per accession, with raw responses cached on
disk so a re-run costs nothing:
  UniProt REST          protein name, gene symbol, sequence length
  AlphaFold DB          latest model version + per-residue pLDDT
  InterPro REST         Pfam domains and their residue ranges
  EBI Proteins API      ClinVar/UniProt clinical significance, protein coordinates

The NMD-escape region is not public annotation: it comes from the project's own
gene table, converted from CDS nucleotides to residues as ceil(nt / 3).
"""
import os, json, math, time, urllib.request, urllib.error

CACHE = "data/cache"
OUT = "site/protein_data"
UA = {"Accept": "application/json"}


def _cache_ok(path):
    """A cached body counts only if it is complete, parseable JSON."""
    try:
        b = open(path, "rb").read()
    except OSError:
        return None
    if b == b"__404__" or b == b"":
        return b                 # InterPro answers "no match" with HTTP 204, empty body
    try:
        json.loads(b)
    except Exception:
        os.remove(path)          # truncated by an interrupted write - refetch
        return None
    return b


def fetch(url, cache_key, retries=3, pause=1.5):
    """GET with on-disk cache. Returns bytes, or None on a definitive 404."""
    path = os.path.join(CACHE, cache_key)
    if os.path.exists(path):
        b = _cache_ok(path)
        if b is not None:
            return None if b == b"__404__" else b
    for k in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                b = r.read()
            tmp = path + ".part"          # atomic: never leave a half file in cache
            open(tmp, "wb").write(b)
            os.replace(tmp, path)
            return b
        except urllib.error.HTTPError as e:
            if e.code == 404:
                open(path, "wb").write(b"__404__")
                return None
            if k == retries - 1:
                raise
            time.sleep(pause * (k + 1))
        except Exception:
            if k == retries - 1:
                raise
            time.sleep(pause * (k + 1))


def uniprot(acc):
    b = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.json", f"up_{acc}.json")
    if b is None:
        return None
    j = json.loads(b)
    desc = j.get("proteinDescription", {})
    nm = (desc.get("recommendedName") or {}).get("fullName", {}).get("value")
    if not nm:
        subs = desc.get("submissionNames") or []
        nm = subs[0]["fullName"]["value"] if subs else acc
    return dict(name=nm, length=j["sequence"]["length"],
                gene=(j.get("genes") or [{}])[0].get("geneName", {}).get("value"))


def alphafold(acc):
    b = fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", f"afmeta_{acc}.json")
    if b is None:
        return None
    meta = json.loads(b)
    if not meta:
        return None
    m = meta[0]
    ver = m.get("latestVersion", 4)
    entry = m.get("entryId", f"AF-{acc}-F1")
    c = fetch(f"https://alphafold.ebi.ac.uk/files/{entry}-confidence_v{ver}.json",
              f"afconf_{acc}.json")
    if c is None:
        return None
    cj = json.loads(c)
    return dict(entry=entry, version=ver, cif=m.get("cifUrl"),
                plddt=[round(float(x), 1) for x in cj["confidenceScore"]])


def pfam(acc):
    b = fetch(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/?page_size=100",
              f"pfam_{acc}.json")
    if not b:
        return []                # None = 404, b"" = HTTP 204, both mean no Pfam match
    out = []
    for r in json.loads(b).get("results", []):
        md = r.get("metadata", {})
        for pr in r.get("proteins", []):
            for loc in pr.get("entry_protein_locations", []):
                for fr in loc.get("fragments", []):
                    out.append(dict(acc=md.get("accession"), name=md.get("name"),
                                    start=fr["start"], end=fr["end"]))
    return sorted(out, key=lambda d: d["start"])


def variants(acc):
    """Pathogenic / likely-pathogenic single-residue missense, grouped by position."""
    b = fetch(f"https://www.ebi.ac.uk/proteins/api/variation/{acc}", f"var_{acc}.json")
    if b is None:
        return [], 0
    feats = json.loads(b).get("features", [])
    per = {}
    total = 0
    for f in feats:
        sig = " ".join(c.get("type", "") for c in f.get("clinicalSignificances", []))
        if "athogenic" not in sig:
            continue
        wt, alt = f.get("wildType"), f.get("alternativeSequence")
        if not wt or not alt or len(wt) != 1 or len(alt) != 1 or alt == "*":
            continue
        try:
            b0, e0 = int(f["begin"]), int(f["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if b0 != e0:
            continue
        total += 1
        d = per.setdefault(b0, dict(pos=b0, wt=wt, alts=set()))
        d["alts"].add(alt)
    rows = [dict(pos=d["pos"], wt=d["wt"], alts="".join(sorted(d["alts"])), n=len(d["alts"]))
            for d in per.values()]
    return sorted(rows, key=lambda d: d["pos"]), total


def nt_to_residue(start_nt, end_nt):
    return int(math.ceil(start_nt / 3)), int(math.ceil(end_nt / 3))


def build_one(rec):
    acc = rec["uniprot"]
    up = uniprot(acc)
    if up is None:
        return None, "uniprot_missing"
    af = alphafold(acc)
    doms = pfam(acc)
    vrows, vtotal = variants(acc)
    esc = nt_to_residue(rec["NMD_region_start"], rec["NMD_region_end"])
    length = up["length"]
    plddt = af["plddt"] if af else []
    if plddt and len(plddt) != length:
        length_note = f"AlphaFold model covers {len(plddt)} of {length} residues"
    else:
        length_note = None
    bands = {"ge90": 0, "b70_90": 0, "b50_70": 0, "lt50": 0}
    for p in plddt:
        k = "ge90" if p >= 90 else "b70_90" if p >= 70 else "b50_70" if p >= 50 else "lt50"
        bands[k] += 1
    return dict(
        key=f"{rec['hgnc_symbol']}_{rec['group']}",
        gene=rec["hgnc_symbol"], uniprot=acc, group=rec["group"],
        transcript=rec["ensembl_transcript_id"],
        name=up["name"], length=length,
        af=None if af is None else dict(entry=af["entry"], version=af["version"], cif=af["cif"]),
        plddt=plddt,
        mean_plddt=round(sum(plddt) / len(plddt), 1) if plddt else None,
        bands=bands, band_note=length_note,
        domains=doms,
        variants=vrows, n_variants=vtotal, n_positions=len(vrows),
        escape=dict(start=esc[0], end=esc[1],
                    cds_start=int(rec["NMD_region_start"]), cds_end=int(rec["NMD_region_end"]),
                    cds_length=int(rec["cds_length"])),
        n_variants_in_escape=sum(r["n"] for r in vrows if esc[0] <= r["pos"] <= esc[1]),
    ), None


def prefetch(records, workers=6, log_every=20):
    """Warm the cache concurrently. The variation endpoint dominates wall time:
    it returns a protein's entire variation record with no server-side filter,
    so it is megabytes per gene while the other three sources are kilobytes.
    Everything lands in the same on-disk cache build_one() reads from."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    os.makedirs(CACHE, exist_ok=True)
    jobs = []
    for r in records:
        a = r["uniprot"]
        jobs += [(f"https://www.ebi.ac.uk/proteins/api/variation/{a}", f"var_{a}.json"),
                 (f"https://rest.uniprot.org/uniprotkb/{a}.json", f"up_{a}.json"),
                 (f"https://alphafold.ebi.ac.uk/api/prediction/{a}", f"afmeta_{a}.json"),
                 (f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{a}/?page_size=100",
                  f"pfam_{a}.json")]
    todo = [(u, k) for u, k in jobs if not os.path.exists(os.path.join(CACHE, k))]
    print(f"prefetch: {len(todo)} of {len(jobs)} responses still to download", flush=True)
    done = errs = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch, u, k): k for u, k in todo}
        for f in as_completed(futs):
            done += 1
            try:
                f.result()
            except Exception as e:
                errs += 1
                print(f"  ! {futs[f]}: {type(e).__name__}", flush=True)
            if done % log_every == 0 or done == len(todo):
                print(f"  prefetch {done}/{len(todo)}, {errs} errors", flush=True)
    # the AlphaFold confidence file needs the version from afmeta, so it is a second wave
    conf = []
    for r in records:
        a = r["uniprot"]
        b = _cache_ok(os.path.join(CACHE, f"afmeta_{a}.json"))
        if not b or b == b"__404__":
            continue
        meta = json.loads(b)
        if not meta:
            continue
        m = meta[0]
        entry, ver = m.get("entryId", f"AF-{a}-F1"), m.get("latestVersion", 4)
        key = f"afconf_{a}.json"
        if not os.path.exists(os.path.join(CACHE, key)):
            conf.append((f"https://alphafold.ebi.ac.uk/files/{entry}-confidence_v{ver}.json", key))
    print(f"prefetch: {len(conf)} confidence files", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(as_completed({ex.submit(fetch, u, k) for u, k in conf}))
    print("prefetch done", flush=True)


def main(records, log_every=10):
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    index, failures = [], []
    for i, rec in enumerate(records, 1):
        try:
            payload, err = build_one(rec)
        except Exception as e:
            payload, err = None, f"{type(e).__name__}: {e}"
        if payload is None:
            failures.append(dict(gene=rec["hgnc_symbol"], uniprot=rec["uniprot"], reason=err))
        else:
            json.dump(payload, open(os.path.join(OUT, payload["key"] + ".json"), "w"))
            index.append({k: payload[k] for k in
                          ("key", "gene", "uniprot", "group", "name", "length",
                           "mean_plddt", "n_variants", "n_positions", "n_variants_in_escape")}
                         | {"n_domains": len(payload["domains"]),
                            "escape_start": payload["escape"]["start"],
                            "escape_end": payload["escape"]["end"],
                            "has_structure": payload["af"] is not None})
        if i % log_every == 0 or i == len(records):
            print(f"{i}/{len(records)} done, {len(failures)} failed", flush=True)
    index.sort(key=lambda d: (d["gene"], d["group"]))
    json.dump(index, open(os.path.join(OUT, "_index.json"), "w"))
    json.dump(failures, open("data/viewer_build_failures.json", "w"))
    return index, failures
