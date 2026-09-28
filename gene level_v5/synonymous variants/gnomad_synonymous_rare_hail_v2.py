"""
gnomAD v4.1 exomes: PASS + synonymous (MANE Select) + AF < 1%, using the public Hail Table.

Local-server version: reads the public AWS S3 mirror of gnomAD anonymously
(no Google or AWS account needed) and writes the output to local disk.

Run inside the conda env:
    conda activate hail
    python gnomad_synonymous_rare_hail.py
"""
import os

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
AF_MAX = 0.01
DRIVER_MEMORY = "16g"   # lower this if the server has less free RAM

# Test on a small region first; set to None for the full exome
# (a full run downloads a very large table and can take many hours).
#INTERVALS = ["chr17:7661779-7687550"]   # TP53
INTERVALS = None

HT_PATH = "s3a://gnomad-public-us-east-1/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht"
OUT = os.path.expanduser("~/gnomad_new/gnomad_exomes_v4.1_synonymous_AFlt1pct.tsv.bgz")
TMP_DIR = os.path.expanduser("~/gnomad_new/hail_tmp")

# Driver memory must be set before Hail/Spark starts
os.environ["PYSPARK_SUBMIT_ARGS"] = f"--driver-memory {DRIVER_MEMORY} pyspark-shell"

import hail as hl  # noqa: E402

os.makedirs(TMP_DIR, exist_ok=True)

hl.init(
    default_reference="GRCh38",
    tmp_dir=TMP_DIR,
    spark_conf={
        # S3 connector matching the Hadoop version bundled with Spark 3.5
        "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
        # Anonymous access to the public gnomAD bucket
        "spark.hadoop.fs.s3a.aws.credentials.provider":
            "org.apache.hadoop.fs.s3a.AnonymousAWSCredentialsProvider",
        "spark.hadoop.fs.s3a.endpoint.region": "us-east-1",
    },
)

# ---------------------------------------------------------------------------
# Analysis (same logic as the original script)
# ---------------------------------------------------------------------------
ht = hl.read_table(HT_PATH)

if INTERVALS:
    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(i, reference_genome="GRCh38") for i in INTERVALS]
    )

# freq[0] = adj genotypes, all samples
ht = ht.filter(
    (hl.len(ht.filters) == 0)            # PASS
    & (ht.freq[0].AC > 0)
    & (ht.freq[0].AF < AF_MAX)
)

# Keep MANE Select transcripts whose only consequence is synonymous_variant
syn = ht.vep.transcript_consequences.filter(
    lambda t: hl.is_defined(t.mane_select)
    & (t.consequence_terms == hl.literal(["synonymous_variant"]))
)
ht = ht.annotate(syn=syn)
ht = ht.filter(hl.len(ht.syn) > 0)
ht = ht.explode("syn")

out = ht.select(
    FILTER="PASS",
    AC=ht.freq[0].AC,
    AN=ht.freq[0].AN,
    AF=ht.freq[0].AF,
    nhomalt=ht.freq[0].homozygote_count,
    gene=ht.syn.gene_symbol,
    gene_id=ht.syn.gene_id,
    transcript=ht.syn.transcript_id,
    HGVSc=ht.syn.hgvsc,
    HGVSp=ht.syn.hgvsp,
)
out.export(OUT)
print(f"\nWrote {OUT}")

# Count from the exported local file rather than out.count(), which would
# re-read the whole gnomAD table from S3 a second time.
n = hl.import_table(OUT, force_bgz=True).count()
print("Variants:", n)

hl.stop()
