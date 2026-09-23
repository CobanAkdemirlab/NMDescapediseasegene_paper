#!/usr/bin/env python3
# Clean workflow schematic (SVG) for NMD-escape disease gene identification
svg = r'''<svg xmlns="http://www.w3.org/2000/svg" width="960" height="1000" viewBox="0 0 960 1000" font-family="Helvetica, Arial, sans-serif">
<defs>
  <marker id="arr" markerWidth="9" markerHeight="9" refX="6" refY="3" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L6,3 L0,6 Z" fill="#3b5b8c"/>
  </marker>
  <marker id="arrO" markerWidth="9" markerHeight="9" refX="6" refY="3" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L6,3 L0,6 Z" fill="#c47b1a"/>
  </marker>
  <marker id="arrP" markerWidth="9" markerHeight="9" refX="6" refY="3" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L6,3 L0,6 Z" fill="#7a52a1"/>
  </marker>
</defs>
<style>
  .title{font-size:23px;font-weight:700;fill:#1a1a1a}
  .h{font-size:17px;font-weight:700}
  .b{font-size:13.5px;fill:#333}
  .i{font-size:12.5px;font-style:italic;fill:#555}
  .big{font-size:26px;font-weight:700;fill:#c47b1a}
</style>

<text x="480" y="34" text-anchor="middle" class="title">Autosomal-dominant NMD-escape disease genes</text>

<!-- 1 ClinVar -->
<rect x="300" y="60" width="360" height="58" rx="9" fill="#dde8f7" stroke="#3b5b8c" stroke-width="1.6"/>
<text x="480" y="84" text-anchor="middle" class="h" fill="#22436e">ClinVar VCF</text>
<text x="480" y="104" text-anchor="middle" class="b">4,281,770 variants</text>

<!-- 2 QC -->
<rect x="270" y="146" width="420" height="150" rx="9" fill="#dde8f7" stroke="#3b5b8c" stroke-width="1.6"/>
<text x="480" y="170" text-anchor="middle" class="h" fill="#22436e">Quality filtering</text>
<text x="300" y="194" class="b">&#8226; GRCh38, autosomal (chr 1&#8211;22)</text>
<text x="300" y="215" class="b">&#8226; valid RefSeq ID</text>
<text x="300" y="236" class="b">&#8226; high-confidence review status</text>
<text x="300" y="257" class="b">&#8226; pathogenic / likely pathogenic only</text>
<text x="300" y="281" class="i">from variant-summary annotation</text>

<!-- 3 NMD-escape annotation -->
<rect x="270" y="324" width="420" height="132" rx="9" fill="#eef0f2" stroke="#6b7280" stroke-width="1.6"/>
<text x="480" y="348" text-anchor="middle" class="h" fill="#374151">NMD-escape annotation</text>
<text x="300" y="372" class="b">&#8226; nonsense variants: aenmd escape calls</text>
<text x="300" y="393" class="b">&#8226; frameshift (+1 / &#8722;1): custom caller,</text>
<text x="318" y="412" class="b">first downstream PTC resolved</text>

<!-- side box -->
<rect x="720" y="330" width="210" height="120" rx="9" fill="#fbe6cf" stroke="#c47b1a" stroke-width="1.8"/>
<text x="825" y="372" text-anchor="middle" class="big">16,493</text>
<text x="825" y="396" text-anchor="middle" class="b">NMD-escape P/LP</text>
<text x="825" y="415" text-anchor="middle" class="b">PTC variants</text>
<text x="825" y="435" text-anchor="middle" class="i">5,266 stop-gain</text>

<!-- 4 enrichment -->
<rect x="270" y="484" width="420" height="112" rx="9" fill="#dde8f7" stroke="#3b5b8c" stroke-width="1.6"/>
<text x="480" y="508" text-anchor="middle" class="h" fill="#22436e">Gene-level enrichment</text>
<text x="480" y="532" text-anchor="middle" class="b">Fisher's exact test, BH FDR &lt; 0.2</text>
<text x="480" y="553" text-anchor="middle" class="b">frameshift +1 / &#8722;1 combined (ACAT)</text>
<text x="480" y="577" text-anchor="middle" class="i">stop-gain and frameshift sets</text>

<!-- 5 result -->
<rect x="255" y="624" width="450" height="60" rx="9" fill="#d9ecd6" stroke="#4a8b45" stroke-width="1.8"/>
<text x="480" y="649" text-anchor="middle" class="h" fill="#2f6b2c">NMD-escape disease genes</text>
<text x="480" y="670" text-anchor="middle" class="i">FDR &lt; 0.05</text>

<!-- 6 comparison -->
<rect x="230" y="712" width="500" height="52" rx="9" fill="#dde8f7" stroke="#3b5b8c" stroke-width="1.6"/>
<text x="480" y="744" text-anchor="middle" class="h" fill="#22436e">Feature comparison vs length-matched controls</text>

<!-- 7 gene / variant -->
<rect x="150" y="800" width="300" height="150" rx="9" fill="#efe7f7" stroke="#7a52a1" stroke-width="1.6"/>
<text x="300" y="826" text-anchor="middle" class="h" fill="#5a3a86">Gene level</text>
<text x="300" y="850" text-anchor="middle" class="b">matched on CDS &amp; escape-</text>
<text x="300" y="869" text-anchor="middle" class="b">region length</text>
<text x="300" y="900" text-anchor="middle" class="b" font-weight="700">&#8594; higher PPI connectivity</text>

<rect x="510" y="800" width="300" height="150" rx="9" fill="#efe7f7" stroke="#7a52a1" stroke-width="1.6"/>
<text x="660" y="826" text-anchor="middle" class="h" fill="#5a3a86">Variant level</text>
<text x="660" y="850" text-anchor="middle" class="b">gene-adjusted mixed model</text>
<text x="660" y="871" text-anchor="middle" class="i">flag ~ is_disease + (1 | transcript)</text>
<text x="660" y="900" text-anchor="middle" class="b" font-weight="700">&#8594; odds ratio per variant</text>

<!-- arrows -->
<line x1="480" y1="118" x2="480" y2="144" stroke="#3b5b8c" stroke-width="2.2" marker-end="url(#arr)"/>
<line x1="480" y1="296" x2="480" y2="322" stroke="#3b5b8c" stroke-width="2.2" marker-end="url(#arr)"/>
<line x1="480" y1="456" x2="480" y2="482" stroke="#3b5b8c" stroke-width="2.2" marker-end="url(#arr)"/>
<line x1="480" y1="596" x2="480" y2="622" stroke="#3b5b8c" stroke-width="2.2" marker-end="url(#arr)"/>
<line x1="480" y1="684" x2="480" y2="710" stroke="#3b5b8c" stroke-width="2.2" marker-end="url(#arr)"/>
<line x1="690" y1="390" x2="718" y2="390" stroke="#c47b1a" stroke-width="2.2" marker-end="url(#arrO)"/>
<!-- split to two -->
<path d="M480,764 L480,782 L300,782 L300,798" fill="none" stroke="#7a52a1" stroke-width="2.2" marker-end="url(#arrP)"/>
<path d="M480,764 L480,782 L660,782 L660,798" fill="none" stroke="#7a52a1" stroke-width="2.2" marker-end="url(#arrP)"/>
</svg>'''
open("schematic.svg","w").write(svg)
# render to PNG
import subprocess, sys
try:
    import cairosvg; cairosvg.svg2png(bytestring=svg.encode(), write_to="schematic.png", output_width=1200, scale=1.0)
    print("rendered with cairosvg")
except Exception as e:
    print("cairosvg missing:", e)
