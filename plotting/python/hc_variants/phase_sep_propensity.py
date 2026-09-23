#!/usr/bin/env python3
"""
Step 4: Phase Separation Propensity Prediction — Full-Length vs Truncated

Computes multiple sequence-based phase separation propensity metrics for
both full-length and truncated (nonsense) versions of each protein, then
compares them to assess whether truncation alters LLPS propensity.

Metrics computed:
  1. Prion-like composition score (PLAAC-style)
     - Fraction of Q/N/G/S/Y residues in sliding windows
  2. Pi-pi interaction propensity (PScore-style)
     - Based on aromatic (Y/F/W) and arginine content driving cation-pi
  3. Charge patterning (SCD — sequence charge decoration)
     - Measures blockiness of charge distribution; higher SCD = more
       charge-separated = more LLPS-prone
  4. Sticker density
     - Fraction of "sticker" residues (Y/F/W/R) that drive multivalent
       LLPS interactions, vs "spacer" residues
  5. Low-complexity fraction
     - Shannon entropy per window; low entropy = more LC = more LLPS-prone

The KEY comparison: does truncation INCREASE phase separation propensity?
If removing a C-terminal IDR solubility tag increases prion-like content,
sticker density, or charge blockiness of the retained protein, that supports
the condensate hypothesis.

Install:
    pip install pandas numpy matplotlib openpyxl requests

Usage:
    python3 phase_sep_propensity.py \
        --hc-variants high_confidence_variants.xlsx \
        --disorder-results ./disorder_results/disorder_analysis_results.csv \
        --charge-results ./charge_results/charge_analysis_results.csv \
        --output-dir ./phase_sep_propensity/

Author: Jacob Schmidt
Date: March 2026
"""

import argparse
import os
import sys
import math
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import requests


# ============================================================================
# Sequence retrieval
# ============================================================================

def fetch_uniprot_sequence(gene_name):
    """Fetch canonical human protein sequence from UniProt by gene name."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'(gene:{gene_name}) AND (organism_id:9606) AND (reviewed:true)',
        'format': 'json', 'size': '1',
        'fields': 'accession,sequence,length',
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code != 200:
            return None, None
        data = resp.json()
        results = data.get('results', [])
        if not results:
            return None, None
        entry = results[0]
        acc = entry.get('primaryAccession', '')
        seq = entry.get('sequence', {}).get('value', '')
        return acc, seq
    except Exception as e:
        print(f"    UniProt error for {gene_name}: {e}")
        return None, None


# ============================================================================
# Phase separation propensity metrics
# ============================================================================

def prion_like_score(seq):
    """
    PLAAC-style prion-like composition score.
    Fraction of residues that are Q, N, G, S, or Y (the canonical
    prion-like domain residues). Higher = more prion-like.
    """
    if not seq:
        return 0.0
    prion_aa = set('QNGSY')
    count = sum(1 for aa in seq if aa in prion_aa)
    return count / len(seq)


def prion_like_score_windowed(seq, window=100):
    """Max prion-like score in any window of given size."""
    if not seq or len(seq) < window:
        return prion_like_score(seq)
    scores = []
    prion_aa = set('QNGSY')
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        scores.append(sum(1 for aa in w if aa in prion_aa) / window)
    return max(scores) if scores else 0.0


def pi_pi_score(seq):
    """
    PScore-style pi-pi / cation-pi interaction propensity.
    Based on aromatic residues (Y, F, W) and arginine (R) content.
    Cation-pi interactions (R-Y, R-F) are the dominant driving force
    for many phase-separating proteins (e.g., FUS, DDX4).
    Score = (fraction aromatic) * (fraction R) * 100
    Higher = more LLPS-prone via pi-pi/cation-pi.
    """
    if not seq:
        return 0.0
    aromatic = sum(1 for aa in seq if aa in 'YFW')
    arg = sum(1 for aa in seq if aa == 'R')
    f_aro = aromatic / len(seq)
    f_arg = arg / len(seq)
    return f_aro * f_arg * 100


def sticker_fraction(seq):
    """
    Fraction of "sticker" residues that drive LLPS multivalent interactions.
    Stickers: Y, F, W (aromatic), R (cation-pi), plus Q and N (polar).
    Based on the sticker-spacer framework (Harmon et al., 2017).
    """
    if not seq:
        return 0.0
    stickers = set('YFWRQN')
    count = sum(1 for aa in seq if aa in stickers)
    return count / len(seq)


def sequence_charge_decoration(seq):
    """
    SCD (Sequence Charge Decoration) metric.
    Measures the blockiness of charge distribution along the sequence.
    Higher SCD = more charge-separated = more LLPS-prone.
    
    SCD = (1/N) * sum_{i<j} (q_i * q_j * sqrt(|j-i|))
    where q_i = +1 for K/R, -1 for D/E, 0 otherwise.
    
    Normalized by N^2 for length-independence.
    """
    if not seq or len(seq) < 2:
        return 0.0
    
    charges = []
    for aa in seq:
        if aa in ('K', 'R'):
            charges.append(1.0)
        elif aa in ('D', 'E'):
            charges.append(-1.0)
        else:
            charges.append(0.0)
    
    N = len(seq)
    scd = 0.0
    # For efficiency, only compute for charged residues
    charged_idx = [(i, charges[i]) for i in range(N) if charges[i] != 0]
    
    for a in range(len(charged_idx)):
        for b in range(a + 1, len(charged_idx)):
            i, qi = charged_idx[a]
            j, qj = charged_idx[b]
            scd += qi * qj * math.sqrt(abs(j - i))
    
    return scd / (N * N)  # normalize by N^2


def shannon_entropy(seq):
    """Per-residue Shannon entropy. Lower = more low-complexity."""
    if not seq:
        return 0.0
    from collections import Counter
    counts = Counter(seq)
    N = len(seq)
    H = 0.0
    for aa, c in counts.items():
        p = c / N
        if p > 0:
            H -= p * math.log2(p)
    return H


def low_complexity_fraction(seq, window=40, entropy_threshold=2.5):
    """
    Fraction of the sequence that falls in low-complexity windows.
    Low complexity = Shannon entropy < threshold in a sliding window.
    """
    if not seq or len(seq) < window:
        return 0.0 if shannon_entropy(seq) >= entropy_threshold else 1.0
    
    lc_residues = [False] * len(seq)
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        if shannon_entropy(w) < entropy_threshold:
            for j in range(i, i + window):
                lc_residues[j] = True
    
    return sum(lc_residues) / len(seq)


def hydrophobic_cluster_score(seq, window=10, threshold=0.45):
    """
    Fraction of sequence in hydrophobic clusters.
    Kyte-Doolittle scale, looking for windows above threshold.
    """
    kd = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }
    if not seq or len(seq) < window:
        return 0.0
    
    in_cluster = [False] * len(seq)
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        avg = np.mean([kd.get(aa, 0) for aa in w])
        if avg > threshold:
            for j in range(i, i + window):
                in_cluster[j] = True
    return sum(in_cluster) / len(seq)


def compute_all_metrics(seq, label=""):
    """Compute all phase separation propensity metrics for a sequence."""
    if not seq:
        return {}
    return {
        'length': len(seq),
        'prion_like_global': prion_like_score(seq),
        'prion_like_max_window': prion_like_score_windowed(seq, 100),
        'pi_pi_score': pi_pi_score(seq),
        'sticker_fraction': sticker_fraction(seq),
        'SCD': sequence_charge_decoration(seq),
        'shannon_entropy': shannon_entropy(seq),
        'low_complexity_frac': low_complexity_fraction(seq),
        'hydrophobic_cluster_frac': hydrophobic_cluster_score(seq),
    }


# ============================================================================
# Main analysis
# ============================================================================

def load_variants(hc_file):
    if hc_file.endswith('.xlsx'):
        df = pd.read_excel(hc_file)
    else:
        df = pd.read_csv(hc_file)
    df = df[df['canonical_is_nonsense'] == True].copy()
    variants = []
    for _, row in df.iterrows():
        variants.append({
            'gene': row['gene'],
            'variant_id': row['variant_id'],
            'aa_pos': int(row['aa_position_canonical']),
            'prot_len': int(row['protein_length_canonical']),
            'position_fraction': row['position_fraction'],
            'pooled_weighted_efficiency': row.get('pooled_weighted_efficiency', np.nan),
        })
    return variants


def main():
    parser = argparse.ArgumentParser(
        description='Phase separation propensity: full-length vs truncated')
    parser.add_argument('--hc-variants', required=True)
    parser.add_argument('--disorder-results', default=None)
    parser.add_argument('--charge-results', default=None)
    parser.add_argument('--output-dir', default='./phase_sep_propensity/')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("Step 4: Phase Separation Propensity — Full-Length vs Truncated")
    print("=" * 70)

    variants = load_variants(args.hc_variants)
    print(f"\nLoaded {len(variants)} variants")

    # Load prior results for cross-referencing
    disorder_df = None
    charge_df = None
    if args.disorder_results and os.path.exists(args.disorder_results):
        disorder_df = pd.read_csv(args.disorder_results)
    if args.charge_results and os.path.exists(args.charge_results):
        charge_df = pd.read_csv(args.charge_results)

    import time
    seen_genes = set()
    all_results = []

    for variant in variants:
        gene = variant['gene']
        if gene in seen_genes:
            continue
        seen_genes.add(gene)

        aa_pos = variant['aa_pos']
        prot_len = variant['prot_len']
        pos_frac = variant['position_fraction']

        print(f"\n--- {gene} (truncation at aa {aa_pos}/{prot_len}, "
              f"{pos_frac:.1%}) ---")

        acc, seq = fetch_uniprot_sequence(gene)
        time.sleep(0.3)
        if not seq:
            print(f"  WARNING: No sequence found, skipping")
            continue
        print(f"  UniProt: {acc}, length: {len(seq)}")

        # Truncated = first aa_pos-1 residues (stop codon at aa_pos)
        trunc_seq = seq[:aa_pos - 1]
        removed_seq = seq[aa_pos - 1:]

        if len(trunc_seq) == 0:
            print(f"  WARNING: Truncation too early, skipping")
            continue

        # Compute metrics
        full_metrics = compute_all_metrics(seq, "full")
        trunc_metrics = compute_all_metrics(trunc_seq, "truncated")
        removed_metrics = compute_all_metrics(removed_seq, "removed")

        # Compute deltas (truncated - full): positive = truncation INCREASES
        deltas = {}
        for key in full_metrics:
            if key == 'length':
                continue
            deltas[f'delta_{key}'] = trunc_metrics[key] - full_metrics[key]

        # Cross-reference
        removed_disorder = np.nan
        pI_shift = np.nan
        if disorder_df is not None:
            match = disorder_df[disorder_df['gene'] == gene]
            if len(match) > 0:
                removed_disorder = match.iloc[0].get('removed_mean_disorder', np.nan)
        if charge_df is not None:
            match = charge_df[charge_df['gene'] == gene]
            if len(match) > 0:
                pI_shift = match.iloc[0].get('pI_shift', np.nan)

        # Flag: does truncation increase LLPS propensity?
        flags = []
        if deltas.get('delta_prion_like_global', 0) > 0.01:
            flags.append('prion_up')
        if deltas.get('delta_sticker_fraction', 0) > 0.01:
            flags.append('sticker_up')
        if deltas.get('delta_SCD', 0) > 0:
            flags.append('SCD_up')
        if deltas.get('delta_pi_pi_score', 0) > 0:
            flags.append('pi_pi_up')
        if deltas.get('delta_low_complexity_frac', 0) < -0.02:
            flags.append('LC_down')  # less LC = less disordered tail buffering
        if deltas.get('delta_hydrophobic_cluster_frac', 0) > 0.01:
            flags.append('hydrophobic_up')

        result = {
            'gene': gene,
            'uniprot_accession': acc,
            'aa_position': aa_pos,
            'protein_length': prot_len,
            'position_fraction': pos_frac,
            'pooled_weighted_efficiency': variant['pooled_weighted_efficiency'],
        }

        # Add all metrics with prefixes
        for key, val in full_metrics.items():
            result[f'full_{key}'] = val
        for key, val in trunc_metrics.items():
            result[f'trunc_{key}'] = val
        for key, val in removed_metrics.items():
            result[f'removed_{key}'] = val
        result.update(deltas)

        result['removed_mean_disorder'] = removed_disorder
        result['pI_shift'] = pI_shift
        result['propensity_flags'] = '; '.join(flags) if flags else 'none'
        result['n_flags'] = len(flags)

        all_results.append(result)

        # Print summary
        print(f"  Prion-like:  full={full_metrics['prion_like_global']:.3f}  "
              f"trunc={trunc_metrics['prion_like_global']:.3f}  "
              f"delta={deltas['delta_prion_like_global']:+.3f}")
        print(f"  Pi-pi:       full={full_metrics['pi_pi_score']:.4f}  "
              f"trunc={trunc_metrics['pi_pi_score']:.4f}  "
              f"delta={deltas['delta_pi_pi_score']:+.4f}")
        print(f"  Sticker:     full={full_metrics['sticker_fraction']:.3f}  "
              f"trunc={trunc_metrics['sticker_fraction']:.3f}  "
              f"delta={deltas['delta_sticker_fraction']:+.3f}")
        print(f"  SCD:         full={full_metrics['SCD']:.6f}  "
              f"trunc={trunc_metrics['SCD']:.6f}  "
              f"delta={deltas['delta_SCD']:+.6f}")
        print(f"  LC frac:     full={full_metrics['low_complexity_frac']:.3f}  "
              f"trunc={trunc_metrics['low_complexity_frac']:.3f}  "
              f"delta={deltas['delta_low_complexity_frac']:+.3f}")
        if flags:
            print(f"  FLAGS: {', '.join(flags)}")

        # ---- Per-gene plot: sliding window comparison ----
        fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
        window = 50
        positions = np.arange(len(seq))

        # 1. Prion-like score
        prion_aa = set('QNGSY')
        prion_profile = []
        for i in range(len(seq)):
            start = max(0, i - window // 2)
            end = min(len(seq), i + window // 2)
            w = seq[start:end]
            prion_profile.append(sum(1 for aa in w if aa in prion_aa) / len(w))

        axes[0].plot(positions, prion_profile, 'b-', alpha=0.7, linewidth=0.8)
        axes[0].axvline(x=aa_pos - 1, color='red', linestyle='--', linewidth=2,
                       label=f'Truncation (aa {aa_pos})')
        axes[0].fill_between(positions, 0, prion_profile,
                            where=positions >= aa_pos - 1,
                            alpha=0.2, color='red', label='Removed')
        axes[0].set_ylabel('Prion-like score\n(Q/N/G/S/Y fraction)')
        axes[0].set_title(f'{gene} — Phase Separation Propensity Profile')
        axes[0].legend(loc='upper right')
        axes[0].set_ylim(0, 1)

        # 2. Sticker density
        sticker_set = set('YFWRQN')
        sticker_profile = []
        for i in range(len(seq)):
            start = max(0, i - window // 2)
            end = min(len(seq), i + window // 2)
            w = seq[start:end]
            sticker_profile.append(sum(1 for aa in w if aa in sticker_set) / len(w))

        axes[1].plot(positions, sticker_profile, 'g-', alpha=0.7, linewidth=0.8)
        axes[1].axvline(x=aa_pos - 1, color='red', linestyle='--', linewidth=2)
        axes[1].fill_between(positions, 0, sticker_profile,
                            where=positions >= aa_pos - 1,
                            alpha=0.2, color='red')
        axes[1].set_ylabel('Sticker fraction\n(Y/F/W/R/Q/N)')
        axes[1].set_ylim(0, 1)

        # 3. Charge (sliding window NCPR)
        charge_map = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
        charge_profile = []
        for i in range(len(seq)):
            start = max(0, i - window // 2)
            end = min(len(seq), i + window // 2)
            w = seq[start:end]
            charge_profile.append(
                sum(charge_map.get(aa, 0) for aa in w) / len(w))

        axes[2].plot(positions, charge_profile, 'purple', alpha=0.7, linewidth=0.8)
        axes[2].axhline(y=0, color='gray', linestyle='-', alpha=0.3)
        axes[2].axvline(x=aa_pos - 1, color='red', linestyle='--', linewidth=2)
        axes[2].fill_between(positions, 0, charge_profile,
                            where=np.array(positions) >= aa_pos - 1,
                            alpha=0.2, color='red')
        axes[2].set_ylabel('NCPR\n(net charge/residue)')
        axes[2].set_xlabel('Amino acid position')
        axes[2].set_ylim(-0.3, 0.3)

        plt.tight_layout()
        plot_file = os.path.join(args.output_dir, f'{gene}_propensity_profile.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()

    # ========================================================================
    # Save results
    # ========================================================================
    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values('n_flags', ascending=False)

    out_csv = os.path.join(args.output_dir, 'phase_sep_propensity_results.csv')
    results_df.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    # ========================================================================
    # Summary comparison plot: delta metrics across all genes
    # ========================================================================
    if len(results_df) > 1:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        
        metrics_to_plot = [
            ('delta_prion_like_global', 'Δ Prion-like score', axes[0, 0]),
            ('delta_pi_pi_score', 'Δ Pi-pi score', axes[0, 1]),
            ('delta_sticker_fraction', 'Δ Sticker fraction', axes[0, 2]),
            ('delta_SCD', 'Δ SCD (charge decoration)', axes[1, 0]),
            ('delta_low_complexity_frac', 'Δ Low-complexity fraction', axes[1, 1]),
            ('delta_hydrophobic_cluster_frac', 'Δ Hydrophobic cluster frac', axes[1, 2]),
        ]

        for col, title, ax in metrics_to_plot:
            if col not in results_df.columns:
                continue
            sorted_df = results_df.sort_values(col, ascending=False)
            colors = ['red' if v > 0 else 'blue' for v in sorted_df[col]]
            ax.barh(range(len(sorted_df)), sorted_df[col], color=colors, alpha=0.7)
            ax.set_yticks(range(len(sorted_df)))
            ax.set_yticklabels(sorted_df['gene'], fontsize=8)
            ax.axvline(x=0, color='black', linewidth=0.5)
            ax.set_title(title, fontsize=10)
            ax.invert_yaxis()

        plt.suptitle('Effect of Truncation on Phase Separation Propensity Metrics\n'
                     '(Red = truncation INCREASES propensity; Blue = decreases)',
                     fontsize=12, fontweight='bold')
        plt.tight_layout()
        summary_plot = os.path.join(args.output_dir, 'delta_propensity_summary.png')
        plt.savefig(summary_plot, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {summary_plot}")

    # ========================================================================
    # Report
    # ========================================================================
    print("\n" + "=" * 70)
    print("PHASE SEPARATION PROPENSITY SUMMARY")
    print("=" * 70)
    print(f"\n{'Gene':12s} {'PosFrac':>8s} {'Δ Prion':>8s} {'Δ Pi-Pi':>8s} "
          f"{'Δ Stick':>8s} {'Δ SCD':>10s} {'Flags':>6s} Notes")
    print("-" * 90)

    for _, row in results_df.iterrows():
        notes = []
        if not np.isnan(row.get('removed_mean_disorder', np.nan)):
            if row['removed_mean_disorder'] > 0.5:
                notes.append(f"rmvd_IDR={row['removed_mean_disorder']:.2f}")
        if not np.isnan(row.get('pI_shift', np.nan)):
            if abs(row['pI_shift']) > 0.3:
                notes.append(f"pI_shift={row['pI_shift']:+.1f}")
        notes_str = '; '.join(notes)

        print(f"{row['gene']:12s} {row['position_fraction']:8.1%} "
              f"{row.get('delta_prion_like_global', 0):+8.3f} "
              f"{row.get('delta_pi_pi_score', 0):+8.4f} "
              f"{row.get('delta_sticker_fraction', 0):+8.3f} "
              f"{row.get('delta_SCD', 0):+10.6f} "
              f"{row['n_flags']:6d} {notes_str}")

    # Interpretation guide
    n_any_increase = (results_df['n_flags'] > 0).sum()
    n_multi = (results_df['n_flags'] >= 2).sum()

    print(f"\nSummary:")
    print(f"  Variants with ANY propensity increase:    {n_any_increase} / {len(results_df)}")
    print(f"  Variants with ≥2 propensity increases:    {n_multi} / {len(results_df)}")

    print(f"\nInterpretation:")
    print(f"  Positive Δ = truncation INCREASES that metric in the retained protein")
    print(f"  Key signals for condensate hypothesis:")
    print(f"    - Δ Prion-like > 0: retained protein is more prion-like (lost non-PrLD tail)")
    print(f"    - Δ SCD > 0: more charge-separated = more LLPS-prone")
    print(f"    - Δ Sticker > 0: higher density of LLPS-driving residues")
    print(f"    - Combined with IDR removal (Step 1) and pI shift (Step 2)")

    print(f"\n{'='*70}")
    print(f"Output files in: {args.output_dir}")
    print(f"  phase_sep_propensity_results.csv  - Full numerical results")
    print(f"  *_propensity_profile.png          - Per-gene sliding window profiles")
    print(f"  delta_propensity_summary.png      - Cross-gene comparison")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()