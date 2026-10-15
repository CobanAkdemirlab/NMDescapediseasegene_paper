"""Replace the viewer's variant track with gene set 4's own disease and control lists.

Source: the variant table Figure 5 is fitted on, one row per variant allele
(`Variant_Key` = chrom:pos|ref|alt, unique across the table). Each row carries
`is_disease`, `type` (fs | snv), the mutation's CDS position and the resulting
PTC position in both nucleotides and amino acids.

  disease  = is_disease 1  -> ClinVar P/LP PTC variants in the case gene
  control  = is_disease 0  -> gnomAD NMD-escape PTC variants in the SAME gene

The contrast is within-gene, so both lists belong on one residue axis: disease
above the backbone, control mirrored below.

Anchor: variants are placed at the MUTATION location, ceil(cds_mutation_loc / 3),
because that is the position Figure 5's variant-level models use. The resulting
PTC residue travels in the tooltip instead - for a stopgain the two coincide, for
a frameshift the PTC lies downstream of the mutation.
"""
import csv, json, math, os, collections

VARIANTS = ("/Users/qkelly/Desktop/track/20261005_set4_bundle_by_figure/"
            "Figure5/input/raresyn_top20_variants.csv")


def load_rows(path=VARIANTS):
    """One entry per variant allele, carrying BOTH anchors.

      m  mutation residue, ceil(cds_mutation_loc / 3)
      p  PTC residue (ptc_aa), or None when the shifted frame yielded no stop
      em / ep  whether that anchor falls in the NMD-escape region

    Membership is evaluated in CDS nucleotides, the native units of both the
    variant position and the region; deciding it on the rounded residue pulls in
    variants that sit 1-2 nt outside but share the boundary codon.

    Both anchors are kept because they answer different questions: Figure 5's
    variant-level models are anchored on the mutation, while the escape region is
    defined by where the PTC lands. For a stopgain the two sit in the same codon;
    for a frameshift the PTC is downstream.
    """
    per = collections.defaultdict(lambda: {"disease": [], "control": []})
    for r in csv.DictReader(open(path, newline="")):
        key = (r["hgnc_symbol"], r["type"])
        arm = "disease" if r["is_disease"] == "1" else "control"
        mut_nt = int(float(r["cds_mutation_loc"]))
        ptc_nt = float(r["ptc_pos"]) if r["ptc_pos"].strip() else None
        ptc_aa = float(r["ptc_aa"]) if r.get("ptc_aa", "").strip() else None
        per[key][arm].append(dict(
            m=int(math.ceil(mut_nt / 3)),
            p=int(ptc_aa) if ptc_aa is not None else None,
            mnt=mut_nt,
            pnt=int(ptc_nt) if ptc_nt is not None else None))
    return per


def load(path=VARIANTS):
    rows = list(csv.DictReader(open(path, newline="")))
    per = collections.defaultdict(lambda: {"disease": {}, "control": {}})
    skipped = []
    for r in rows:
        key = (r["hgnc_symbol"], r["type"])
        arm = "disease" if r["is_disease"] == "1" else "control"
        try:
            loc = int(float(r["cds_mutation_loc"]))
        except (ValueError, KeyError):
            skipped.append((r.get("Variant_Key"), "no cds_mutation_loc"))
            continue
        res = int(math.ceil(loc / 3))
        ptc = r.get("ptc_aa", "").strip()
        ptc_res = int(float(ptc)) if ptc else None
        slot = per[key][arm].setdefault(res, {"pos": res, "n": 0, "ptc": [], "no_ptc": 0,
                                              "nts": []})
        slot["n"] += 1
        slot["nts"].append(loc)
        if ptc_res is None:
            slot["no_ptc"] += 1
        elif ptc_res not in slot["ptc"]:
            slot["ptc"].append(ptc_res)
    out = {}
    for key, arms in per.items():
        out[key] = {arm: sorted(v.values(), key=lambda d: d["pos"]) for arm, v in arms.items()}
    return out, len(rows), skipped


def attach2(site_dir, per=None):
    """Write per-variant rows with both anchors, plus totals for each."""
    per = per if per is not None else load_rows()
    pdir = os.path.join(site_dir, "protein_data")
    index = json.load(open(os.path.join(pdir, "_index.json")))
    tot = collections.Counter()
    for row in index:
        path = os.path.join(pdir, row["key"] + ".json")
        d = json.load(open(path))
        ns, ne = d["escape"]["cds_start"], d["escape"]["cds_end"]
        rec = {}
        for arm in ("disease", "control"):
            out = []
            for v in per.get((d["gene"], d["group"]), {}).get(arm, []):
                out.append(dict(m=v["m"], p=v["p"],
                                em=bool(ns <= v["mnt"] <= ne),
                                ep=None if v["pnt"] is None else bool(ns <= v["pnt"] <= ne)))
            rec[arm] = out
            rec["n_" + arm] = len(out)
            rec["n_%s_in_escape_mutation" % arm] = sum(1 for v in out if v["em"])
            rec["n_%s_in_escape_ptc" % arm] = sum(1 for v in out if v["ep"])
            rec["n_%s_no_ptc" % arm] = sum(1 for v in out if v["p"] is None)
            tot[arm] += len(out)
            tot[arm + "_esc_m"] += rec["n_%s_in_escape_mutation" % arm]
            tot[arm + "_esc_p"] += rec["n_%s_in_escape_ptc" % arm]
        d["nmdesc_variants"] = rec
        json.dump(d, open(path, "w"))
        for k in ("n_disease", "n_control"):
            row[k] = rec[k]
        row["n_disease_in_escape"] = rec["n_disease_in_escape_mutation"]
        row["n_control_in_escape"] = rec["n_control_in_escape_mutation"]
    json.dump(index, open(os.path.join(pdir, "_index.json"), "w"))
    return tot


def attach(site_dir, table=None):
    table = table or load()[0]
    pdir = os.path.join(site_dir, "protein_data")
    index = json.load(open(os.path.join(pdir, "_index.json")))
    stats = collections.Counter()
    for row in index:
        path = os.path.join(pdir, row["key"] + ".json")
        d = json.load(open(path))
        arms = table.get((d["gene"], d["group"]), {"disease": [], "control": []})
        # membership is decided in CDS nucleotides - the native units of both the
        # variant position and the region. Deciding it on the rounded residue
        # pulls in variants that sit 1-2 nt outside but share the boundary codon.
        ns, ne = d["escape"]["cds_start"], d["escape"]["cds_end"]
        rec = {}
        for arm in ("disease", "control"):
            lst = []
            n_in = 0
            for x in arms.get(arm, []):
                inside = sum(1 for nt in x["nts"] if ns <= nt <= ne)
                n_in += inside
                lst.append(dict(pos=x["pos"], n=x["n"], ptc=sorted(x["ptc"]),
                                no_ptc=x["no_ptc"], n_in_escape=inside))
            rec[arm] = lst
            rec["n_" + arm] = sum(x["n"] for x in lst)
            rec["n_%s_pos" % arm] = len(lst)
            rec["n_%s_in_escape" % arm] = n_in
        d["nmdesc_variants"] = rec
        json.dump(d, open(path, "w"))
        row["n_disease"] = rec["n_disease"]
        row["n_control"] = rec["n_control"]
        row["n_disease_in_escape"] = rec["n_disease_in_escape"]
        row["n_control_in_escape"] = rec["n_control_in_escape"]
        stats["disease"] += rec["n_disease"]
        stats["control"] += rec["n_control"]
        stats["records"] += 1
    json.dump(index, open(os.path.join(pdir, "_index.json"), "w"))
    return stats


if __name__ == "__main__":
    t = attach2("restore/NMDesc-protein-viewer")
    print(f"disease {t['disease']}  in escape: by mutation {t['disease_esc_m']}, "
          f"by PTC {t['disease_esc_p']}")
    print(f"control {t['control']}  in escape: by mutation {t['control_esc_m']}, "
          f"by PTC {t['control_esc_p']}")
    print(f"total {t['disease'] + t['control']}")
