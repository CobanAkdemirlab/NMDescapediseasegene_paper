"""ParSe v2, ported line-for-line from stevewhitten/Parse_v2_FASTA js/parse.js
+ js/part-b.js (the 'PS potential' tab) at js/constants.js defaults
WINDOW_SIZE=25, REGION_LENGTH=20, CUTOFF=0.9.

Returns the two quantities the NMDesc Figure 6E model consumes:
  p_sum          'Σ classifier distance P'   = Σ dist_norm over region=='P' residues
  longest_ps_idr 'length of longest PS IDR'  = longest P region, reported only if >=20

Faithful to the shipped code, including two quirks:
  * U_pi / U_q are inert in the shipped build (the pi-pi and SCD block is commented
    out, so scd=0 and U_pi=0); only `region`/`dist_norm` are used here, not the
    `_pi_q` variants, matching the plain 'Σ classifier distance P' column.
  * the terminal floor(w/2) residues inherit `region` from the first/last scored
    window but keep dist_norm = 0 (the copy lines are commented out upstream),
    so they extend regions without adding to p_sum.
"""
import math

AA = {
    'A': (0.37, 1.42, 0.0728), 'C': (0.25, 0.73, 0.3557), 'D': (0.30, 1.01, -0.0552),
    'E': (0.42, 1.63, -0.0295), 'F': (0.17, 1.16, 0.4201), 'G': (0.13, 0.50, -0.0589),
    'H': (0.20, 1.20, 0.0874), 'I': (0.39, 1.12, 0.3805), 'K': (0.56, 1.24, -0.0053),
    'L': (0.24, 1.29, 0.3819), 'M': (0.36, 1.21, 0.1613), 'N': (0.27, 0.71, -0.0390),
    'P': (1.00, 0.65, -0.0492), 'Q': (0.53, 1.02, 0.0126), 'R': (0.38, 1.06, 0.0394),
    'S': (0.24, 0.71, -0.0282), 'T': (0.32, 0.78, 0.0239), 'V': (0.39, 0.99, 0.2947),
    'W': (0.25, 1.05, 0.4114), 'Y': (0.25, 0.67, 0.3113), 'X': (0.37, 1.42, 0.0728),
    'U': (0.25, 0.73, 0.3557),
}
WINDOW_SIZE, REGION_LENGTH, CUTOFF = 25, 20, 0.9
C1, C2 = -0.244078945, 0.7885823
_M = -1.0 / C1


def _anchor(c3, c4):
    b = c3 - _M * c4
    x = (b - C2) / (C1 - _M)
    y = _M * x + b
    return math.sqrt((c4 - x) ** 2 + (c3 - y) ** 2)


ID_DIST = _anchor(0.5582901, 1.022552)
PS_DIST = _anchor(0.5416, 0.9327272)


def residues(seq, w=WINDOW_SIZE):
    n = len(seq)
    region = [None] * n
    dist = [0.0] * n
    half = w // 2
    for i in range(n - w + 1):
        mid = i + half
        win = seq[i:i + w]
        ppii = helix = hydr = 0.0
        nD = nE = nK = nR = 0
        for ch in win:
            p = AA.get(ch)
            if p is None:
                continue
            ppii += p[0]; helix += p[1]; hydr += p[2]
            if ch == 'D': nD += 1
            elif ch == 'E': nE += 1
            elif ch == 'K': nK += 1
            elif ch == 'R': nR += 1
        ppii /= w; helix /= w; hydr /= w
        if hydr >= 0.08280152:
            region[mid] = 'F'
            dist[mid] = (hydr - 0.08280152) / (0.01679414 * 2.0)
            continue
        net_charge = abs(nD + nE - nK - nR)
        if ppii == 1.0:
            ppii = 0.98
        v_exp = 0.503 - 0.11 * math.log(1.0 - ppii)
        rh = 2.16 * (4 * w) ** v_exp + 0.26 * 4 * net_charge - 0.29 * (4 * w) ** 0.5
        nu = math.log(rh / 2.16) / math.log(4 * w)
        b = nu - _M * helix
        x = (b - C2) / (C1 - _M)
        y = _M * x + b
        d = math.sqrt((helix - x) ** 2 + (nu - y) ** 2)
        if ((nu - C2) / C1) <= helix:
            region[mid] = 'D'; dist[mid] = d / ID_DIST
        else:
            region[mid] = 'P'; dist[mid] = d / PS_DIST
    for i in range(half):                       # terminal inheritance: region only
        region[i] = region[half]
        region[n - i - 1] = region[n - half - 1]
    return region, dist


def _regions(region, label, cutoff=CUTOFF):
    out, i, n = [], 0, len(region)
    while i < n:
        j, found, cp, cw = i, False, 0, 0
        while j < n:
            cw += 1
            if region[j] == label:
                cp += 1
            if cp / cw >= cutoff:
                found = True
            else:
                break
            j += 1
        if found:
            out.append((i, j - 1)); i = j
        else:
            i += 1
    return out


def parse_v2(seq):
    """-> (p_sum, longest_ps_idr) or (None, None) when ParSe refuses the sequence."""
    s = ''.join(str(seq).split()).upper()
    if len(s) < 25 or len(s) > 10000:
        return None, None
    region, dist = residues(s)
    p_sum = round(sum(d for r, d in zip(region, dist) if r == 'P'), 1)
    ps = _regions(region, 'P')
    longest = max((e - st for st, e in ps), default=0)
    return p_sum, (longest if longest >= REGION_LENGTH else 0)
