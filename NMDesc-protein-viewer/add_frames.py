"""Split the frameshift design by reading frame: per-frame escape regions, and a
frame label on every frameshift variant.

Two sources:
  PTC_info20260201_region.csv   one row per (transcript, plus1|plus2) with that
                                frame's own can_region_start / can_region_end
  the gene set 4 variant table  one row per variant allele, keyed chrom:pos|ref|alt

FRAME OF A VARIANT. Deleting one base shifts everything downstream LEFT by one,
so the downstream sequence is read in the +1 offset of the original coordinates;
inserting one base shifts it right, giving +2. The pipeline's frame therefore
follows (len(ref) - len(alt)) mod 3, not the net length change:

    net -1 nt (deletion of 1)  -> pipeline 'plus1'  -> reported here as "-1"
    net +1 nt (insertion of 1) -> pipeline 'plus2'  -> reported here as "+1"
    net  0 mod 3               -> in-frame indel, not a frameshift, label None

That direction is not assumed. It was established against the per-frame regions:
restricted to the variants where the two frame regions do NOT overlap - the only
informative cases - the variant's own PTC falls in the region of the frame this
mapping assigns 360 times against 91 for the opposite assignment.

The two frame regions are also checked against the gene table: round(median) of
the two reproduces NMD_region_start / NMD_region_end on all 94 frameshift
records, which is how build_gene_all.R derives it.
"""
import csv, json, math, os, collections

REGIONS = "/Users/qkelly/Desktop/clinvar/PTC_info20260201_region.csv"
VARIANTS = ("/Users/qkelly/Desktop/track/20261005_set4_bundle_by_figure/"
            "Figure5/input/raresyn_top20_variants.csv")
# pipeline frame -> the label a geneticist uses for it
FRAME_LABEL = {"plus1": "-1", "plus2": "+1"}


def frame_of(variant_key):
    """-> 'plus1' | 'plus2' | None (in-frame indel or unparseable key)."""
    parts = variant_key.split("|")
    if len(parts) != 3:
        return None
    _, ref, alt = parts
    return {1: "plus1", 2: "plus2"}.get((len(ref) - len(alt)) % 3)


def load_regions(path=REGIONS):
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(path, newline="")):
        try:
            s, e = int(float(r["can_region_start"])), int(float(r["can_region_end"]))
        except (ValueError, TypeError, KeyError):
            continue
        out[r["transcript"]][r["type"]] = {"cds_start": s, "cds_end": e,
                                           "start": int(math.ceil(s / 3)),
                                           "end": int(math.ceil(e / 3))}
    return out


def load_variant_frames(path=VARIANTS):
    """(gene, 'fs') -> {Variant_Key: 'plus1'|'plus2'|None}"""
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(path, newline="")):
        if r["type"] != "fs":
            continue
        out[(r["hgnc_symbol"], "fs")][r["Variant_Key"]] = frame_of(r["Variant_Key"])
    return out


def key_order(path=VARIANTS):
    """Variant keys in the same order attach2 wrote them, per (gene, type, arm)."""
    out = collections.defaultdict(list)
    for r in csv.DictReader(open(path, newline="")):
        arm = "disease" if r["is_disease"] == "1" else "control"
        out[(r["hgnc_symbol"], r["type"], arm)].append(r["Variant_Key"])
    return out


def attach(site_dir):
    regions, vframes, korder = load_regions(), load_variant_frames(), key_order()
    pdir = os.path.join(site_dir, "protein_data")
    index = json.load(open(os.path.join(pdir, "_index.json")))
    stats = collections.Counter()
    for row in index:
        path = os.path.join(pdir, row["key"] + ".json")
        d = json.load(open(path))
        if d["group"] == "fs":
            fr = regions.get(d["transcript"], {})
            if "plus1" in fr and "plus2" in fr:
                d["escape_frames"] = {FRAME_LABEL[k]: fr[k] for k in ("plus1", "plus2")}
                stats["records_with_frames"] += 1
            else:
                d["escape_frames"] = None
                stats["records_without_frames"] += 1
            NV = d.get("nmdesc_variants")
            if NV:
                for arm in ("disease", "control"):
                    keys = korder.get((d["gene"], "fs", arm), [])
                    if len(keys) != len(NV[arm]):
                        raise SystemExit(
                            f"{row['key']}/{arm}: {len(NV[arm])} variants but "
                            f"{len(keys)} keys - ordering assumption broken")
                    for v, vk in zip(NV[arm], keys):
                        f = vframes[(d["gene"], "fs")].get(vk)
                        v["f"] = FRAME_LABEL.get(f)        # None for in-frame indels
                        stats["frame_" + str(v["f"])] += 1
        else:
            d["escape_frames"] = None
            NV = d.get("nmdesc_variants")
            if NV:
                for arm in ("disease", "control"):
                    for v in NV[arm]:
                        v["f"] = None
        json.dump(d, open(path, "w"))
    json.dump(index, open(os.path.join(pdir, "_index.json"), "w"))
    return stats


if __name__ == "__main__":
    print(dict(attach("restore/NMDesc-protein-viewer")))
