# NMDesc protein viewer — gene set 4 case genes

A static, browser-only viewer for the 158 case records of NMDesc gene set 4
(`4_raresyn_top20`): 94 frameshift + 64 stopgain, 134 distinct genes (24 genes
appear in both designs with different NMD-escape regions, so each design gets
its own record).

Four panels per protein: the AlphaFold model coloured by per-residue pLDDT, the
pLDDT distribution, a Pfam domain track with a pathogenic-missense lollipop, and
the positions carrying the most distinct substitutions. The **NMD-escape region**
is drawn on the track and marked in the table — that layer is the project's own
definition and exists in no public tool.

## Run it

`fetch()` is blocked on `file://`, so serve the folder over HTTP:

```bash
cd site
python3 -m http.server 8000
# open http://localhost:8000/
```

To publish on GitHub Pages, commit this folder as `docs/` on the default branch
and enable Pages → "main / docs".

## Layout

```
site/
  index.html                  the whole application, no build step, no framework
  protein_data/
    _index.json               158 rows, drives the gene list and the search box
    <GENE>_<fs|snv>.json      one self-contained payload per record
```

`index.html` loads one payload per selection. Nothing else is bundled: the 3D
coordinates stream from AlphaFold DB at view time, which keeps the whole data
directory at 1.3 MB.

## Payload fields

| field | meaning |
|---|---|
| `gene`, `uniprot`, `group`, `transcript` | identity; `group` is `fs` or `snv` |
| `name`, `length` | UniProt recommended name, residue count |
| `af` | AlphaFold entry id, model version, mmCIF URL (`null` when AFDB has no model) |
| `plddt` | per-residue pLDDT, one value per residue |
| `mean_plddt`, `bands` | mean, and residue counts in the four confidence bands |
| `domains` | Pfam signature matches: accession, name, start, end |
| `nmdesc_variants` | gene set 4's own lists: `disease[]` / `control[]`, each entry `{pos, n, ptc[], no_ptc, n_in_escape}`, plus per-arm totals |

| `escape` | NMD-escape region in CDS nucleotides and in residues |
| `n_variants_in_escape` | how many of `n_variants` fall inside that region |
| `tracks` | Table S4 feature lanes: `domains_s4`, `features`, `slims`, `morfs`, `ptms`, `nls`, `lcs`; each a list of `{label, start, end}` |
| `exons` | exon structure in residues, the last exon-exon junction, and whether it is safe to draw |
| `s4_protein` | the Table S4 protein name this record was matched to, or `null` |
| `flags` | the ten gene-level flags for this gene in this design, `true` / `false` / `null` |

## Feature lanes and gene-level flags

The lanes under the variant track come from **Banani et al., Supplemental Table S4**
(Dev Cell 2022;57(14):1776-1788.e8, doi:10.1016/j.devcel.2022.06.010, PMID 35809564,
PMC9339523) — the same source the project's own flags are
computed from, so the lanes and the flags cannot disagree about what was annotated.
Positions were parsed out of the table's `name (start-end, start-end, ...)` cells.

The flag rule is worth stating because it is easy to misread. `annotate_motif_flags.R`
computes

```
flag = max(position of that feature class) * 3 >= NMD_region_start
```

so a flag is TRUE when the **most C-terminal** feature of the class *reaches* the start
of the escape region. It is **not** a test of overlap with the region, and it ignores
every feature except the last one. A record can therefore carry a TRUE flag with no
feature block drawn inside the teal band.

`NA` means Table S4 has no annotation of that class for the protein. Two classes are
annotated so sparsely across this gene set that their flags carry almost no signal:

| flag | TRUE | FALSE | NA |
|---|---:|---:|---:|
| `SLiM_flag` | 7 | 8 | 143 |
| `MoRF_flag` | 1 | 1 | 156 |

`NLS_flag` is also mostly absent (31 / 9 / 118). The remaining flags are populated for
nearly every record.

## Exon structure and the 50-nucleotide rule

Exons come from Ensembl (`/lookup/id/{transcript}?expand=1`) for the transcript named
in the project gene table. Only the coding part of each exon is used; it is placed on
the residue axis as `ceil(CDS nt / 3)`, so two exons share a residue when the boundary
falls inside a codon. The dashed marker is the **last exon-exon junction**.

That junction is the point the NMD-escape definition is anchored to, and the data agree
with the definition exactly:

| design | records | escape start == last junction − 50 nt |
|---|---:|---|
| stopgain (`snv`) | 64 | **64 / 64, zero offset** |
| frameshift (`fs`) | 94 | 1 / 94 — median offset −75 nt |

The stopgain column reproduces the 50-nucleotide rule on every record, which is an
independent check on both the exon mapping here and the escape coordinates in the gene
table. Frameshift records are expected to differ: that design uses the plus1/plus2
region definition, not the 50-nt rule.

**Eleven records do not get an exon track.** For these the Ensembl transcript named in the gene
table and the UniProt accession beside it are different isoforms, so the protein lengths disagree
by more than one residue. Drawing exon boundaries from one isoform onto the other's residue axis
would put every boundary in the wrong place, so the lane states the mismatch instead. The other
panels for those records are unaffected.

| gene | design | transcript | UniProt | Ensembl aa | UniProt aa |
|---|---|---|---|---:|---:|
| BPTF | fs | `ENST00000306378` | Q12830 | 2921 | 3046 |
| EZH2 | fs | `ENST00000320356` | Q15910 | 752 | 746 |
| FGA | fs | `ENST00000403106` | P02671 | 645 | 866 |
| GATA4 | fs | `ENST00000532059` | P43694 | 444 | 442 |
| HNRNPK | fs | `ENST00000376263` | P61978 | 465 | 463 |
| NKX2-1 | fs | `ENST00000354822` | P43699 | 402 | 371 |
| RARB | fs | `ENST00000330688` | P10826 | 449 | 455 |
| RARB | snv | `ENST00000330688` | P10826 | 449 | 455 |
| RORB | fs | `ENST00000376896` | Q92753 | 460 | 470 |
| THRA | fs | `ENST00000450525` | P10827 | 411 | 490 |
| TRPS1 | fs | `ENST00000395715` | Q9UHF7 | 1295 | 1281 |

This is a pointer at the gene table, not only at the viewer: wherever these pairings feed a
length- or position-dependent step, the transcript and the accession are describing different
proteins. `exon_isoform_mismatches.csv` in the run folder carries the same eleven rows.

## Sources

| source | endpoint | gives |
|---|---|---|
| UniProt | `rest.uniprot.org/uniprotkb/{acc}.json` | name, sequence length |
| AlphaFold DB | `alphafold.ebi.ac.uk/api/prediction/{acc}` + `-confidence_v{n}.json` | model version, per-residue pLDDT |
| InterPro | `ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/` | Pfam domains and ranges |
| gene set 4 variant table | `Figure5/input/raresyn_top20_variants.csv` | disease and control PTC variants, CDS and residue coordinates |

Variants kept are single-residue missense classified **Pathogenic** or **Likely
pathogenic**. The size of a lollipop is the number of *distinct substitutions*
at that position, not the number of submissions.

ClinVar was read through the EBI Proteins API rather than NCBI E-utilities:
E-utilities requires a contact email that was not available, and the EBI endpoint
carries the same ClinVar classifications already mapped to protein coordinates,
which is the form the track needs.

The NMD-escape region comes from the project gene table, converted from CDS
nucleotides to residues as `ceil(nt / 3)`.

## Coverage and known gaps

- **158 / 158** records built, 0 failures.
- **155 / 158** have an AlphaFold model. BPTF (3046 aa), DYNC1H1 (4646 aa) and
  PKD1 (4303 aa) return 404 from AlphaFold DB — they exceed the length limit for
  a single-chain model. Their other three panels still work.
- **3 records** have no Pfam match (GATM, MN1 in both designs). InterPro answers
  HTTP 204 for these; the domain track renders the bare backbone.
- **3 records** have no pathogenic missense in ClinVar/UniProt (HCK, NOVA2,
  ODC1). The lollipop track is empty for them, which is a real absence, not a
  fetch failure.
- 14,966 pathogenic missense variants across the set, 3,309 of them inside an
  NMD-escape region. That ratio is a raw count over heterogeneous genes, not a
  tested enrichment — no statistical claim is made here or in the interface.

## Rebuilding

`build_viewer_data.py` regenerates `protein_data/` from the gene table. Every raw
HTTP response is cached under `data/cache/`, so a re-run with a warm cache costs
no network; `prefetch(records, workers=6)` warms it concurrently.

## A note on the Table S4 citation

Some copies of this data circulate with the citation `Cell. 2023;184(2):341-359.e25, doi:10.1016/j.cell.2023.01.013`. That is wrong. The title, authors and PMCID in it are right, but the journal, year, volume, pages and DOI are not, and that DOI resolves to an unrelated paper (*A molecular device for the redox quality control of GroEL/ES substrates*, Cell 2023;186:1039-1049.e17). Checked against Europe PMC (`PMCID:PMC9339523`) and Crossref (`10.1016/j.devcel.2022.06.010`); both return Developmental Cell 2022;57(14):1776-1788.e8.

## Why the escape band ends before the last exon (frameshift records)

The two designs use different intervals, and `build_gene_all.R` is explicit about it.

**Stopgain (`snv`)** — a transcript interval:

```r
nmdesc_end   = max(cds_end)                                    # end of the CDS
nmdesc_start = max(cds_end) - 50 - last_exon_length            # 50 nt before the last exon
```

so the band always runs to the end of the CDS: 64 / 64 records have `NMD_region_end == cds_length`, and the start reproduces the 50-nt rule exactly on all 64.

**Frameshift (`fs`)** — not a transcript interval at all:

```r
median_can_region_start = round(median(can_region_start))   # over the +1 / +2 alleles
median_can_region_end   = round(median(can_region_end))     # of that transcript
```

Each frameshift allele has its own PTC-to-new-stop region, and the gene-level region is the
**median across those alleles**. Its right edge is the median position of the *frame-shifted* stop
codon, which sits upstream of the wild-type stop — so the band ends before the last exon does, by
construction. None of the 94 frameshift records reaches the CDS end (median shortfall 82.5 nt,
range 11-659), and because the edge is a rounded median over variants rather than a structural
boundary it need not be codon-aligned: 73 / 94 are not.

Practical consequence: the teal band means *the last coding exon plus 50 nt* on a stopgain record
and *where a typical frameshift allele stops translating* on a frameshift record. They are not the
same quantity and should not be compared edge to edge.

## Variants: gene set 4's own disease and control lists

The track shows the table Figure 5 is fitted on, not generic public annotation. One row per
variant allele (`Variant_Key` = `chrom:pos|ref|alt`, unique across all 2,541 rows).

| arm | definition | n | fs | snv |
|---|---|---:|---:|---:|
| disease | ClinVar P/LP PTC variants in the case gene | 1731 | 1318 | 413 |
| control | gnomAD NMD-escape PTC variants in the SAME gene | 810 | 595 | 215 |

The contrast is **within-gene** - the control arm is other variants in the same protein - so there
is no matching ratio here, unlike the gene-level figures. Disease variants are drawn above the
backbone, control mirrored below, on one shared scale.

**Anchor is selectable.** A lollipop sits at either the MUTATION location (`ceil(cds_mutation_loc / 3)`, the anchor Figure 5's variant-level models use) or the PTC location (`ptc_aa`, where translation actually stops). The dropdown above the track switches between them; the other coordinate always travels in the tooltip and in the table's second column. Each payload stores one row per variant allele with both anchors and both escape flags, so switching is a regroup, not a refetch.

For a stopgain the two coincide to within a codon — the mutated base is one of the three bases of the new stop (offsets 0, 1 or 2 nt over 625 of 628 stopgain rows). For a frameshift the PTC lies downstream, sometimes far downstream. The choice therefore changes the escape tally:

| arm | n | inside escape, by mutation | by PTC |
|---|---:|---:|---:|
| disease | 1731 | 1514 (87%) | 1475 (86%) |
| control | 810 | 699 (86%) | 660 (83%) |

21 frameshift rows (8 disease, 13 control) have no PTC recorded; they appear on the mutation view and are omitted — with a count — from the PTC view.

Three stopgain rows, all NKX2-5 `snv_disease`, carry `ptc_pos` 973 = `cds_end` − 2, i.e. the NATURAL stop rather than a premature one, while their mutations sit at CDS 414, 433 and 451. The PTC was evidently not resolved for them and fell back to the CDS end; on the PTC view all three stack on the last codon. Worth checking in the source table.

**Several alleles can share one position.** 138 mutation positions carry more than one allele
(up to 7, e.g. MSH6 fs at CDS 3953, where seven different indels at the same base give seven
different PTCs). The lollipop size is the number of alleles, so these are not duplicates.

**Escape membership is decided in CDS nucleotides**, the native units of both the variant position
and the region - not on the rounded residue, which pulls in 8 variants that sit 1-2 nt outside but
share the boundary codon. By that rule 1514 / 1731 disease (87%) and 699 / 810 control (86%) fall
inside. The remainder fall outside because for frameshift records the gene-level region is a median
over alleles, so individual alleles land on either side of it.

**Gaps.** 39 records have no control variant (no gnomAD NMD-escape PTC in that gene) and the lower
arm is empty. One record, SLC40A1 `fs`, has no disease variant in this table.

The earlier ClinVar/UniProt missense set from the EBI Proteins API has been removed from the
payloads. It was filtered with a substring test on `clinicalSignificances` that would also have
matched `Conflicting interpretations of pathogenicity`; that was never audited, so it is not
shipped rather than shipped unverified.

## Reading frames in the frameshift design

The frameshift escape region is the median of two frames, so both are drawn: a **+1** band, a **−1** band, and the median band the gene-level analysis actually uses. Per-frame regions come from `PTC_info20260201_region.csv`, one row per (transcript, plus1|plus2).

**Check:** `round(median(plus1, plus2))` reproduces `NMD_region_start` / `NMD_region_end` on **94 / 94** frameshift records, which is exactly how `build_gene_all.R` derives the gene-level region. The two tables are consistent and the join is right.

**Frame of a variant.** Deleting one base shifts the downstream sequence LEFT, so it is read in the +1 offset of the original coordinates; inserting one base gives +2. The pipeline frame therefore follows `(len(ref) - len(alt)) mod 3`, which is the OPPOSITE sign to the net length change:

| label here | net length change | pipeline frame | n |
|---|---|---|---:|
| `+1` | +1 nt (insertion) | `plus2` | 850 |
| `−1` | −1 nt (deletion) | `plus1` | 1052 |
| none | multiple of 3 | — in-frame indel | 11 |

That direction was not assumed. Restricted to variants where the two frame regions do NOT overlap — the only cases that discriminate — the variant's own PTC falls inside the region of the frame this mapping assigns **360** times, against **91** for the opposite assignment. Over all frameshift variants, 89% land inside their assigned frame's region.

Attempting to validate instead against the `PTC_loc` column failed: its positions do not line up with `ptc_pos` under any constant offset, so the two columns use different coordinate conventions. The region-based check above is what the frame labels rest on.

**Marker shape encodes the frame** so colour stays free for the disease / control arm: circle = +1, diamond = −1, square = in-frame indel. Variants at one residue in different frames are kept as separate markers — they stop at different places and must not be merged.
