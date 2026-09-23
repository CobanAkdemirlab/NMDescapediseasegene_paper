#!/usr/bin/env python3
"""
Step 2: Charge Distribution Analysis for NMD-Escaping Truncation Variants

Builds on the disorder analysis (Step 1) to examine whether removing
C-terminal IDRs shifts the protein's charge profile in ways that could
promote aberrant biomolecular condensate formation.

Specifically tests:
  1. Is the removed region enriched in acidic (D/E) residues vs the retained protein?
     (Loss of acidic IDR = loss of "solubility tag" → promotes condensation)
  2. Does truncation shift the protein's estimated pI toward basic?
     (Basic shift → risk of nucleolar mislocalization, as per Mensah et al. HMGB1)
  3. Does the charge patterning (blocks of same-sign charge) change?
     (Charge blocks drive electrostatic phase separation)
  4. What is the NCPR (net charge per residue) and FCR (fraction charged residues)
     in the Das-Pappu phase diagram framework?

Inputs:
    --hc-variants       high_confidence_variants.xlsx
    --disorder-results  disorder_analysis_results.csv from Step 1 (optional,
                        used to cross-reference disorder scores)
    --output-dir        Output directory

Install:
    pip install pandas matplotlib numpy requests openpyxl

Usage:
    python3 charge_analysis.py \
        --hc-variants high_confidence_variants.xlsx \
        --output-dir ./charge_results/

Author: Jacob Schmidt
Date: March 2026
"""

import argparse
import os
import sys
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import requests


# ============================================================================
# Amino acid properties
# ============================================================================

# Charge at pH 7
CHARGE_PH7 = {'R': 1, 'K': 1, 'H': 0.1, 'D': -1, 'E': -1}

# pKa values for pI estimation (N-term, C-term, and charged sidechains)
PKA = {
    'Nterm': 9.69, 'Cterm': 2.34,
    'D': 3.65, 'E': 4.25, 'C': 8.18,
    'Y': 10.07, 'H': 6.00, 'K': 10.54, 'R': 12.48,
}

# Kyte-Doolittle
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'X': 0.0, 'U': 0.0,
}


# ============================================================================
# Isoelectric point estimation
# ============================================================================

def charge_at_pH(sequence, pH):
    """Calculate net charge of a protein at a given pH."""
    charge = 0.0
    # N-terminus (positive)
    charge += 1.0 / (1.0 + 10 ** (pH - PKA['Nterm']))
    # C-terminus (negative)
    charge += -1.0 / (1.0 + 10 ** (PKA['Cterm'] - pH))

    for aa in sequence:
        if aa in ('D', 'E', 'C', 'Y'):
            pka = PKA[aa]
            charge += -1.0 / (1.0 + 10 ** (pka - pH))
        elif aa in ('H', 'K', 'R'):
            pka = PKA[aa]
            charge += 1.0 / (1.0 + 10 ** (pH - pka))

    return charge


def estimate_pI(sequence, precision=0.01):
    """Estimate isoelectric point by bisection."""
    low, high = 0.0, 14.0
    while (high - low) > precision:
        mid = (low + high) / 2.0
        c = charge_at_pH(sequence, mid)
        if c > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


# ============================================================================
# Das-Pappu phase diagram parameters
# ============================================================================

def compute_das_pappu(sequence):
    """
    Compute the fraction of charged residues (FCR) and net charge per
    residue (NCPR) used in the Das-Pappu phase separation framework.

    Returns (FCR, NCPR, f+, f-)
    """
    n = len(sequence)
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan

    n_pos = sum(1 for aa in sequence if aa in ('R', 'K'))
    n_neg = sum(1 for aa in sequence if aa in ('D', 'E'))

    f_pos = n_pos / n
    f_neg = n_neg / n
    FCR = f_pos + f_neg
    NCPR = f_pos - f_neg

    return FCR, NCPR, f_pos, f_neg


# ============================================================================
# Charge block analysis
# ============================================================================

def find_charge_blocks(sequence, min_block_len=5, min_density=0.5):
    """
    Identify contiguous blocks of predominantly positive or negative charge.

    A block is defined as a window of at least min_block_len residues where
    the fraction of same-sign charged residues >= min_density.

    Returns list of (start, end, sign, density) tuples.
    """
    n = len(sequence)
    blocks = []

    for sign, target_aas in [('positive', ('R', 'K')), ('negative', ('D', 'E'))]:
        # Sliding window
        for wsize in range(min_block_len, min(n + 1, 51)):
            for start in range(n - wsize + 1):
                window = sequence[start:start + wsize]
                count = sum(1 for aa in window if aa in target_aas)
                density = count / wsize
                if density >= min_density:
                    # Check not already covered by a larger block
                    is_new = True
                    for existing in blocks:
                        if (existing[2] == sign and
                                existing[0] <= start and existing[1] >= start + wsize):
                            is_new = False
                            break
                    if is_new:
                        # Remove smaller blocks fully contained in this one
                        blocks = [b for b in blocks if not (
                            b[2] == sign and b[0] >= start and b[1] <= start + wsize)]
                        blocks.append((start, start + wsize, sign, density))

    # Sort by position
    blocks.sort(key=lambda x: x[0])
    return blocks


# ============================================================================
# Sequence fetching (same as step 1)
# ============================================================================

def lookup_uniprot_id(gene_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'(gene:{gene_name}) AND (organism_id:9606) AND (reviewed:true)',
        'format': 'json', 'size': '5',
        'fields': 'accession,gene_names,length',
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code != 200:
            return None, None
        data = resp.json()
        results = data.get('results', [])
        if not results:
            return None, None
        return results[0]['primaryAccession'], results[0].get('sequence', {}).get('length', 0)
    except Exception as e:
        print(f"    UniProt lookup error for {gene_name}: {e}")
        return None, None


def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            lines = resp.text.strip().split('\n')
            return ''.join(lines[1:])
        return None
    except:
        return None


def load_variants(hc_file, min_pos_frac=0.0):
    if hc_file.endswith('.xlsx'):
        df = pd.read_excel(hc_file)
    else:
        df = pd.read_csv(hc_file)
    df = df[df['canonical_is_nonsense'] == True].copy()
    if min_pos_frac > 0:
        df = df[df['position_fraction'] >= min_pos_frac]
    variants = []
    for _, row in df.iterrows():
        variants.append({
            'gene': row['gene'],
            'variant_id': row['variant_id'],
            'variant_type': row.get('variant_type', 'unknown'),
            'aa_pos': int(row['aa_position_canonical']),
            'prot_len': int(row['protein_length_canonical']),
            'position_fraction': row['position_fraction'],
            'pooled_weighted_efficiency': row.get('pooled_weighted_efficiency', np.nan),
        })
    return variants


# ============================================================================
# Plotting
# ============================================================================

def plot_charge_profile(gene, sequence, aa_pos, output_dir):
    """Plot sliding-window charge profile with truncation point."""
    n = len(sequence)
    positions = np.arange(1, n + 1)
    window = 21

    # Per-residue charge
    charges = np.array([CHARGE_PH7.get(aa, 0) for aa in sequence])
    charge_smooth = np.convolve(charges, np.ones(window) / window, mode='same')

    # Cumulative charge from C-terminus (shows what is lost by truncation)
    charges_from_c = np.cumsum(charges[::-1])[::-1]

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    # Panel 1: Sliding window charge
    ax = axes[0]
    ax.fill_between(positions, charge_smooth, where=charge_smooth >= 0,
                    alpha=0.3, color='blue')
    ax.fill_between(positions, charge_smooth, where=charge_smooth < 0,
                    alpha=0.3, color='red')
    ax.plot(positions, charge_smooth, color='gray', linewidth=0.5)
    ax.axvline(x=aa_pos, color='darkred', linewidth=2, label=f'Truncation (aa {aa_pos})')
    ax.axvspan(aa_pos, n, alpha=0.15, color='red')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel(f'Charge Density\n(window={window})')
    ax.set_title(f'{gene} \u2014 Charge Distribution', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9)

    # Panel 2: Cumulative charge from C-terminus
    ax = axes[1]
    ax.plot(positions, charges_from_c, color='purple', linewidth=1)
    ax.axvline(x=aa_pos, color='darkred', linewidth=2)
    ax.axvspan(aa_pos, n, alpha=0.15, color='red')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel('Cumulative Charge\n(from C-terminus)')
    ax.set_xlabel('Amino Acid Position')

    # Annotate the total charge lost
    charge_lost = sum(CHARGE_PH7.get(aa, 0) for aa in sequence[aa_pos:])
    ax.annotate(f'Charge removed: {charge_lost:+.1f}',
                xy=(aa_pos, charges_from_c[aa_pos] if aa_pos < n else 0),
                fontsize=10, color='darkred',
                xytext=(aa_pos - n * 0.15, charges_from_c[min(aa_pos, n - 1)] + 3),
                arrowprops=dict(arrowstyle='->', color='darkred'))

    plt.tight_layout()
    outpath = os.path.join(output_dir, f'{gene}_charge_profile.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {outpath}")


def plot_pI_shift(results_df, output_dir):
    """Bar chart showing pI of full-length vs truncated protein."""
    fig, ax = plt.subplots(figsize=(14, 6))
    genes = results_df['gene'].values
    x = np.arange(len(genes))

    ax.bar(x - 0.2, results_df['pI_full'], 0.35, label='Full-length protein',
           color='steelblue', alpha=0.8)
    ax.bar(x + 0.2, results_df['pI_truncated'], 0.35, label='Truncated protein',
           color='indianred', alpha=0.8)

    # Draw arrows showing direction of shift
    for i, (_, row) in enumerate(results_df.iterrows()):
        shift = row['pI_shift']
        if abs(shift) > 0.05:
            color = 'red' if shift > 0 else 'blue'
            ax.annotate(f'{shift:+.2f}', (i, max(row['pI_full'], row['pI_truncated']) + 0.15),
                        ha='center', fontsize=8, color=color, fontweight='bold')

    ax.axhline(y=7.0, color='gray', linestyle='--', alpha=0.5, label='Neutral pH')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Estimated Isoelectric Point (pI)')
    ax.set_title('Effect of Truncation on Protein Isoelectric Point',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'pI_shift_comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_das_pappu(results_df, output_dir):
    """
    Plot full-length vs truncated proteins on the Das-Pappu diagram
    (FCR vs NCPR). Arrows show the shift caused by truncation.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Background regions (simplified Das-Pappu)
    ax.axhline(y=0, color='gray', linewidth=0.5)
    ax.axvline(x=0, color='gray', linewidth=0.5)

    # Plot each variant as full-length -> truncated with arrow
    for _, row in results_df.iterrows():
        # Full length point
        ax.scatter(row['NCPR_full'], row['FCR_full'],
                   color='steelblue', s=60, zorder=5)
        # Truncated point
        ax.scatter(row['NCPR_truncated'], row['FCR_truncated'],
                   color='indianred', s=60, zorder=5)
        # Arrow from full to truncated
        dx = row['NCPR_truncated'] - row['NCPR_full']
        dy = row['FCR_truncated'] - row['FCR_full']
        if abs(dx) > 0.001 or abs(dy) > 0.001:
            ax.annotate('', xy=(row['NCPR_truncated'], row['FCR_truncated']),
                        xytext=(row['NCPR_full'], row['FCR_full']),
                        arrowprops=dict(arrowstyle='->', color='gray', lw=1))
        # Label
        ax.annotate(row['gene'],
                    (row['NCPR_truncated'], row['FCR_truncated']),
                    fontsize=8, ha='left', va='bottom',
                    xytext=(5, 3), textcoords='offset points')

    # Legend
    ax.scatter([], [], color='steelblue', s=60, label='Full-length')
    ax.scatter([], [], color='indianred', s=60, label='Truncated')
    ax.legend(fontsize=10)

    ax.set_xlabel('Net Charge Per Residue (NCPR)', fontsize=12)
    ax.set_ylabel('Fraction Charged Residues (FCR)', fontsize=12)
    ax.set_title('Das-Pappu Charge Space: Full-Length vs Truncated',
                 fontsize=14, fontweight='bold')

    plt.tight_layout()
    outpath = os.path.join(output_dir, 'das_pappu_diagram.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Charge distribution analysis for NMD-escaping truncation variants')
    parser.add_argument('--hc-variants', required=True,
                        help='Path to high_confidence_variants.xlsx')
    parser.add_argument('--disorder-results', default=None,
                        help='Path to disorder_analysis_results.csv from Step 1 (optional)')
    parser.add_argument('--output-dir', default='./charge_results/',
                        help='Output directory')
    parser.add_argument('--min-position-fraction', type=float, default=0.0,
                        help='Only analyze variants with position_fraction >= this value')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("Step 2: Charge Distribution Analysis")
    print("NMD-Escaping Truncation Variants — Condensate Hypothesis")
    print("=" * 70)

    # Load variants
    print(f"\nLoading variants from {args.hc_variants}...")
    variants = load_variants(args.hc_variants, args.min_position_fraction)
    print(f"  {len(variants)} variants loaded")

    # Optionally load disorder results for cross-referencing
    disorder_df = None
    if args.disorder_results and os.path.exists(args.disorder_results):
        disorder_df = pd.read_csv(args.disorder_results)
        print(f"  Loaded disorder results from {args.disorder_results}")

    # Process each variant
    gene_cache = {}
    all_results = []

    for variant in variants:
        gene = variant['gene']
        aa_pos = variant['aa_pos']
        prot_len = variant['prot_len']

        print(f"\n--- {gene} ({variant['variant_id']}) ---")
        print(f"  Truncation: aa {aa_pos}/{prot_len} ({variant['position_fraction']:.1%})")

        # Fetch sequence
        if gene in gene_cache:
            uniprot_id, sequence = gene_cache[gene]
        else:
            uniprot_id, _ = lookup_uniprot_id(gene)
            if uniprot_id is None:
                print(f"  SKIPPING: No UniProt accession found")
                continue
            time.sleep(0.5)
            sequence = fetch_uniprot_sequence(uniprot_id)
            if sequence is None:
                print(f"  SKIPPING: Could not fetch sequence")
                continue
            gene_cache[gene] = (uniprot_id, sequence)
            time.sleep(0.5)

        actual_len = len(sequence)
        if aa_pos >= actual_len:
            print(f"  SKIPPING: Truncation beyond sequence length")
            continue

        full_seq = sequence
        trunc_seq = sequence[:aa_pos]
        removed_seq = sequence[aa_pos:]

        # pI estimation
        pI_full = estimate_pI(full_seq)
        pI_trunc = estimate_pI(trunc_seq)
        pI_shift = pI_trunc - pI_full

        # Das-Pappu parameters
        FCR_full, NCPR_full, fpos_full, fneg_full = compute_das_pappu(full_seq)
        FCR_trunc, NCPR_trunc, fpos_trunc, fneg_trunc = compute_das_pappu(trunc_seq)
        FCR_rm, NCPR_rm, fpos_rm, fneg_rm = compute_das_pappu(removed_seq)

        # Amino acid composition of removed region
        n_rm = len(removed_seq)
        n_D = removed_seq.count('D')
        n_E = removed_seq.count('E')
        n_R = removed_seq.count('R')
        n_K = removed_seq.count('K')
        net_charge_rm = n_R + n_K - n_D - n_E

        # Charge blocks in removed region
        blocks = find_charge_blocks(removed_seq) if n_rm >= 5 else []
        n_neg_blocks = sum(1 for b in blocks if b[2] == 'negative')
        n_pos_blocks = sum(1 for b in blocks if b[2] == 'positive')

        # Cross-reference disorder if available
        removed_disorder = np.nan
        if disorder_df is not None:
            match = disorder_df[disorder_df['variant_id'] == variant['variant_id']]
            if len(match) > 0:
                removed_disorder = match.iloc[0].get('removed_mean_disorder', np.nan)

        result = {
            'gene': gene,
            'variant_id': variant['variant_id'],
            'variant_type': variant['variant_type'],
            'uniprot_accession': uniprot_id,
            'protein_length': actual_len,
            'truncation_position': aa_pos,
            'removed_residues': n_rm,
            'position_fraction': aa_pos / actual_len,
            'pooled_weighted_efficiency': variant['pooled_weighted_efficiency'],
            # pI
            'pI_full': pI_full,
            'pI_truncated': pI_trunc,
            'pI_shift': pI_shift,
            'pI_shift_direction': 'basic' if pI_shift > 0.1 else ('acidic' if pI_shift < -0.1 else 'neutral'),
            # Das-Pappu full
            'FCR_full': FCR_full,
            'NCPR_full': NCPR_full,
            'frac_pos_full': fpos_full,
            'frac_neg_full': fneg_full,
            # Das-Pappu truncated
            'FCR_truncated': FCR_trunc,
            'NCPR_truncated': NCPR_trunc,
            'frac_pos_truncated': fpos_trunc,
            'frac_neg_truncated': fneg_trunc,
            # Das-Pappu removed
            'FCR_removed': FCR_rm,
            'NCPR_removed': NCPR_rm,
            'frac_pos_removed': fpos_rm,
            'frac_neg_removed': fneg_rm,
            # Removed region composition
            'removed_n_Asp': n_D,
            'removed_n_Glu': n_E,
            'removed_n_Arg': n_R,
            'removed_n_Lys': n_K,
            'removed_net_charge': net_charge_rm,
            'removed_charge_per_residue': net_charge_rm / n_rm if n_rm > 0 else np.nan,
            # Charge blocks
            'removed_n_neg_blocks': n_neg_blocks,
            'removed_n_pos_blocks': n_pos_blocks,
            # Disorder cross-ref
            'removed_mean_disorder': removed_disorder,
            # Removed sequence
            'removed_sequence': removed_seq,
        }

        # Condensate risk flags
        flags = []
        if pI_shift > 0.3:
            flags.append(f"pI_SHIFT_BASIC (+{pI_shift:.2f})")
        if n_rm > 0 and (n_D + n_E) / n_rm > 0.15:
            flags.append(f"REMOVED_ACIDIC_RICH ({(n_D+n_E)/n_rm:.0%} D/E)")
        if NCPR_rm < -0.05:
            flags.append("REMOVED_NET_NEGATIVE")
        if not np.isnan(removed_disorder) and removed_disorder > 0.5:
            flags.append("REMOVED_IS_IDR")

        result['condensate_charge_flags'] = '; '.join(flags) if flags else 'none'
        all_results.append(result)

        # Print
        print(f"  UniProt: {uniprot_id}")
        print(f"  pI full: {pI_full:.2f}  pI truncated: {pI_trunc:.2f}  "
              f"shift: {pI_shift:+.2f} ({result['pI_shift_direction']})")
        print(f"  Das-Pappu full:      FCR={FCR_full:.3f}  NCPR={NCPR_full:+.4f}")
        print(f"  Das-Pappu truncated: FCR={FCR_trunc:.3f}  NCPR={NCPR_trunc:+.4f}")
        print(f"  Removed region ({n_rm} aa): D={n_D} E={n_E} R={n_R} K={n_K}  "
              f"net={net_charge_rm:+d}")
        print(f"  Charge flags: {result['condensate_charge_flags']}")

        # Plot
        plot_charge_profile(gene, sequence, aa_pos, args.output_dir)

    # Save results
    if not all_results:
        print("\nERROR: No results. Check connectivity.")
        sys.exit(1)

    results_df = pd.DataFrame(all_results)

    csv_cols = [c for c in results_df.columns if c != 'removed_sequence']
    out_csv = os.path.join(args.output_dir, 'charge_analysis_results.csv')
    results_df[csv_cols].to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    # Summary plots
    print("\nGenerating summary plots...")
    plot_pI_shift(results_df, args.output_dir)
    plot_das_pappu(results_df, args.output_dir)

    # Report
    print("\n" + "=" * 70)
    print("CHARGE ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\n{'Gene':10s} {'Remvd':>6s} {'pI_full':>8s} {'pI_trunc':>8s} "
          f"{'pI_shift':>8s} {'NCPR_rm':>8s} {'D+E_rm':>7s} {'R+K_rm':>7s} {'Flags'}")
    print("-" * 105)
    for _, row in results_df.iterrows():
        n_rm = row['removed_residues']
        print(f"{row['gene']:10s} {n_rm:6d} "
              f"{row['pI_full']:8.2f} {row['pI_truncated']:8.2f} "
              f"{row['pI_shift']:+8.2f} {row['NCPR_removed']:+8.4f} "
              f"{row['removed_n_Asp']+row['removed_n_Glu']:7d} "
              f"{row['removed_n_Arg']+row['removed_n_Lys']:7d} "
              f"{row['condensate_charge_flags']}")

    n_basic = (results_df['pI_shift'] > 0.3).sum()
    n_acidic_rich = results_df['condensate_charge_flags'].str.contains('ACIDIC_RICH').sum()

    print(f"\nCondensate-relevant charge findings:")
    print(f"  Variants with basic pI shift (>0.3):     {n_basic} / {len(results_df)}")
    print(f"  Variants removing acidic-rich region:     {n_acidic_rich} / {len(results_df)}")

    print(f"\n{'='*70}")
    print(f"Output files in: {args.output_dir}")
    print(f"  charge_analysis_results.csv   - Full numerical results")
    print(f"  <GENE>_charge_profile.png     - Per-gene charge profiles")
    print(f"  pI_shift_comparison.png       - pI shift bar chart")
    print(f"  das_pappu_diagram.png         - FCR vs NCPR phase diagram")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()