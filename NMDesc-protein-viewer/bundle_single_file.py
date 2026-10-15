"""Fold the viewer's JSON payloads into index.html so the page opens from file://.

Browsers refuse fetch() on the file:// protocol, so the served version cannot read
protein_data/*.json when the HTML is simply double-clicked. This script embeds the
index and every payload in a <script type="application/json"> block and installs a
shim that answers the app's own fetch() calls from that block. Requests the shim
does not recognise - notably the AlphaFold coordinate download - fall through to
the real fetch untouched.

The app source is not modified: the shim is injected ahead of it.
"""
import json, os, re

SHIM = """
<script id="nmdesc-embedded-data" type="application/json">__DATA__</script>
<script>
/* Serve the app's own JSON requests from the embedded block so the page works
   from file://. Anything else (AlphaFold coordinates) goes to the network. */
(function () {
  var DATA = JSON.parse(document.getElementById('nmdesc-embedded-data').textContent);
  var realFetch = (typeof window.fetch === 'function') ? window.fetch.bind(window) : null;
  function canned(obj) {
    return Promise.resolve({
      ok: true, status: 200,
      json: function () { return Promise.resolve(obj); },
      text: function () { return Promise.resolve(JSON.stringify(obj)); }
    });
  }
  window.fetch = function (url, opts) {
    var m = String(url).match(/protein_data\\/(.+)\\.json$/);
    if (m) {
      var obj = (m[1] === '_index') ? DATA.index : DATA.payloads[m[1]];
      if (obj !== undefined) return canned(obj);
      return Promise.reject(new Error('not embedded: ' + m[1]));
    }
    if (realFetch) return realFetch(url, opts);
    return Promise.reject(new Error('fetch unavailable'));
  };
})();
</script>
"""

OFFLINE_NOTE = """
<div style="margin:0 0 14px;padding:9px 12px;background:#FFF8E6;border:1px solid #F0DFB0;
     border-radius:9px;font-size:11.5px;color:#7A5C14">
  Standalone build: all 158 records are embedded in this file, so it works offline and from
  <code>file://</code>. The 3D panel still streams coordinates from AlphaFold DB and needs
  internet; every other panel is local.
</div>
"""


def build(site_dir, out_path):
    html = open(os.path.join(site_dir, "index.html"), encoding="utf-8").read()
    pdir = os.path.join(site_dir, "protein_data")
    index = json.load(open(os.path.join(pdir, "_index.json"), encoding="utf-8"))
    payloads = {}
    for fn in sorted(os.listdir(pdir)):
        if fn.endswith(".json") and fn != "_index.json":
            payloads[fn[:-5]] = json.load(open(os.path.join(pdir, fn), encoding="utf-8"))
    missing = [r["key"] for r in index if r["key"] not in payloads]
    if missing:
        raise SystemExit(f"index references {len(missing)} payloads not on disk: {missing[:5]}")

    blob = json.dumps({"index": index, "payloads": payloads}, separators=(",", ":"))
    # a literal </script> inside the block would close it early
    blob = blob.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    shim = SHIM.replace("__DATA__", blob)

    if "<script src=" not in html:
        raise SystemExit("could not find the 3Dmol script tag to inject before")
    html = html.replace("<script src=", shim + "<script src=", 1)
    html = html.replace('<main id="main">', '<main id="main">' + OFFLINE_NOTE, 1)
    html = html.replace("<title>NMDesc protein viewer",
                        "<title>NMDesc protein viewer (standalone)", 1)
    open(out_path, "w", encoding="utf-8").write(html)
    return len(index), len(payloads), os.path.getsize(out_path)


if __name__ == "__main__":
    n_ix, n_pl, size = build("restore/NMDesc-protein-viewer",
                             "NMDesc_protein_viewer_standalone.html")
    print(f"{n_ix} index rows, {n_pl} payloads embedded, {size/1e6:.2f} MB")
