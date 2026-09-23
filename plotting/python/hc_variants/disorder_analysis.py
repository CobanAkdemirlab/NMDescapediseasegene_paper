#!/usr/bin/env python3
"""
Intrinsic Disorder Analysis of C-Terminal Truncation Regions

Reads variants directly from high_confidence_variants.xlsx, maps gene names
to UniProt accessions via the UniProt API, fetches protein sequences, and
runs per-residue disorder prediction on each.

Install requirements:
    pip install metapredict pandas matplotlib numpy requests openpyxl

Usage:
    # Fetch sequences only (for PARSEv2 input) — no disorder analysis:
    python3 disorder_analysis.py \
        --hc-variants high_confidence_variants.xlsx \
        --fetch-fasta tier1_uniprot.fasta

    # Analyze all 22 high-confidence variants:
    python3 disorder_analysis.py \
        --hc-variants high_confidence_variants.xlsx \
        --output-dir ./disorder_results/

    # Only analyze near-C-terminal truncations (position_fraction >= 0.85):
    python3 disorder_analysis.py \
        --hc-variants high_confidence_variants.xlsx \
        --output-dir ./disorder_results/ \
        --min-position-fraction 0.85

Inputs:
    --hc-variants    Path to high_confidence_variants.xlsx
    --fetch-fasta    If set, just fetch UniProt sequences to this FASTA
                     file and exit (no disorder analysis). Use the output
                     with make_parsev2_fasta.py.

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
import requests

try:
    import metapredict as meta
    USE_METAPREDICT = True
    print("Using metapredict for disorder prediction")
except ImportError:
    USE_METAPREDICT = False
    print("WARNING: metapredict not available; using amino acid propensity method")


# ============================================================================
# Amino acid property scales
# ============================================================================

HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'X': 0.0, 'U': 0.0, 'B': -3.5, 'Z': -3.5,
}

CHARGE = {'R': 1, 'K': 1, 'H': 0.1, 'D': -1, 'E': -1}

DISORDER_PROPENSITY = {
    'A': 0.06, 'R': 0.18, 'N': 0.13, 'D': 0.19, 'C': -0.20,
    'Q': 0.18, 'E': 0.30, 'G': 0.17, 'H': 0.05, 'I': -0.49,
    'L': -0.34, 'K': 0.26, 'M': -0.10, 'F': -0.42, 'P': 0.41,
    'S': 0.14, 'T': 0.05, 'W': -0.49, 'Y': -0.34, 'V': -0.38,
    'X': 0.0, 'U': 0.0,
}


# ============================================================================
# Load variants
# ============================================================================

def load_variants(hc_file, min_position_fraction=0.0):
    if hc_file.endswith('.xlsx'):
        df = pd.read_excel(hc_file)
    else:
        df = pd.read_csv(hc_file)

    df = df[df['canonical_is_nonsense'] == True].copy()

    if min_position_fraction > 0:
        df = df[df['position_fraction'] >= min_position_fraction]

    variants = []
    for _, row in df.iterrows():
        variants.append({
            'gene': row['gene'],
            'variant_id': row['variant_id'],
            'variant_type': row.get('variant_type', 'unknown'),
            'aa_pos': int(row['aa_position_canonical']),
            'prot_len': int(row['protein_length_canonical']),
            'position_fraction': row['position_fraction'],
            'enst_mane_select': row.get('enst_mane_select', ''),
            'enst_canonical': row.get('enst_canonical', ''),
            'pooled_weighted_efficiency': row.get('pooled_weighted_efficiency', np.nan),
            'tier': row.get('tier', ''),
            'gnomad_pLI': row.get('gnomad_pLI', np.nan),
            'gnomad_LOEUF': row.get('gnomad_LOEUF', np.nan),
        })

    print(f"  Loaded {len(variants)} variants from {hc_file}")
    return variants


# ============================================================================
# UniProt fetching
# ============================================================================

def lookup_uniprot_id(gene_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'(gene:{gene_name}) AND (organism_id:9606) AND (reviewed:true)',
        'format': 'json',
        'size': '5',
        'fields': 'accession,gene_names,length',
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code != 200:
            print(f"    UniProt search returned status {resp.status_code} for {gene_name}")
            return None, None
        data = resp.json()
        results = data.get('results', [])
        if not results:
            print(f"    No UniProt results for {gene_name}")
            return None, None
        entry = results[0]
        accession = entry['primaryAccession']
        length = entry.get('sequence', {}).get('length', 0)
        return accession, length
    except Exception as e:
        print(f"    UniProt lookup error for {gene_name}: {e}")
        return None, None


def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            lines = resp.text.strip().split('\n')
            header = lines[0]
            seq = ''.join(lines[1:])
            return seq, header
        else:
            print(f"    UniProt FASTA returned status {resp.status_code} for {uniprot_id}")
            return None, None
    except Exception as e:
        print(f"    Could not fetch sequence for {uniprot_id}: {e}")
        return None, None


# ============================================================================
# --fetch-fasta mode: just get sequences, skip all disorder analysis
# ============================================================================

def fetch_fasta_only(hc_file, out_fasta, min_position_fraction=0.0):
    """
    Fetch canonical UniProt sequences for all variants in the HC xlsx and
    write them to a multi-entry FASTA file. No disorder analysis is run.
    The output is ready for make_parsev2_fasta.py.
    """
    print("=" * 70)
    print("Fetching UniProt sequences for PARSEv2 input")
    print("=" * 70)

    variants = load_variants(hc_file, min_position_fraction)
    if not variants:
        print("ERROR: No variants loaded.")
        return

    seen_genes = {}
    written = 0

    with open(out_fasta, 'w') as fh:
        for variant in variants:
            gene = variant['gene']

            if gene in seen_genes:
                print(f"  {gene}: already fetched (skipping duplicate)")
                continue

            print(f"\n  {gene}")
            uniprot_id, uniprot_len = lookup_uniprot_id(gene)
            if uniprot_id is None:
                print(f"    SKIPPED: no UniProt accession found")
                continue
            print(f"    UniProt: {uniprot_id}  ({uniprot_len} aa)")
            time.sleep(0.5)

            sequence, header = fetch_uniprot_sequence(uniprot_id)
            if sequence is None:
                print(f"    SKIPPED: sequence fetch failed")
                continue

            if abs(len(sequence) - variant['prot_len']) > 10:
                print(f"    NOTE: fetched length {len(sequence)} vs expected "
                      f"{variant['prot_len']} — check isoform")

            seen_genes[gene] = (uniprot_id, sequence)

            # UniProt-format header so make_parsev2_fasta.py parses accession correctly
            fh.write(f">sp|{uniprot_id}|{gene}_HUMAN {gene} OS=Homo sapiens\n")
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i+60] + "\n")
            written += 1
            time.sleep(0.5)

    print(f"\nWrote {written} sequences to {out_fasta}")
    print(f"\nNext step:")
    print(f"  python make_parsev2_fasta.py --fasta {out_fasta}")


# ============================================================================
# Disorder prediction
# ============================================================================

def predict_disorder_metapredict(sequence):
    scores = meta.predict_disorder(sequence)
    return np.array(scores)


def predict_disorder_propensity(sequence, window=21):
    scores = np.array([DISORDER_PROPENSITY.get(aa, 0.0) for aa in sequence])
    half = window // 2
    smoothed = np.zeros(len(sequence))
    for i in range(len(sequence)):
        start = max(0, i - half)
        end = min(len(sequence), i + half + 1)
        smoothed[i] = np.mean(scores[start:end])
    smoothed = (smoothed - smoothed.min()) / (smoothed.max() - smoothed.min() + 1e-10)
    return smoothed


def predict_disorder(sequence):
    if USE_METAPREDICT:
        return predict_disorder_metapredict(sequence)
    else:
        return predict_disorder_propensity(sequence)


# ============================================================================
# Sequence analysis
# ============================================================================

def analyze_region(sequence):
    if len(sequence) == 0:
        return {k: np.nan for k in [
            'length', 'mean_hydrophobicity', 'net_charge', 'charge_per_residue',
            'frac_charged', 'frac_positive', 'frac_negative',
            'frac_aromatic', 'frac_proline', 'frac_glycine',
        ]}
    n = len(sequence)
    hydro = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in sequence])
    charges = [CHARGE.get(aa, 0) for aa in sequence]
    net_charge = sum(charges)
    n_pos = sum(1 for aa in sequence if aa in ('R', 'K'))
    n_neg = sum(1 for aa in sequence if aa in ('D', 'E'))
    n_aro = sum(1 for aa in sequence if aa in ('F', 'W', 'Y'))
    n_pro = sum(1 for aa in sequence if aa == 'P')
    n_gly = sum(1 for aa in sequence if aa == 'G')
    return {
        'length': n,
        'mean_hydrophobicity': hydro,
        'net_charge': net_charge,
        'charge_per_residue': net_charge / n,
        'frac_charged': (n_pos + n_neg) / n,
        'frac_positive': n_pos / n,
        'frac_negative': n_neg / n,
        'frac_aromatic': n_aro / n,
        'frac_proline': n_pro / n,
        'frac_glycine': n_gly / n,
    }


# ============================================================================
# Plotting
# ============================================================================

def plot_disorder_profile(gene, sequence, disorder_scores, aa_pos, output_dir):
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1, 1]})
    positions = np.arange(1, len(sequence) + 1)
    window = 21

    ax = axes[0]
    ax.fill_between(positions, disorder_scores, alpha=0.3, color='steelblue')
    ax.plot(positions, disorder_scores, color='steelblue', linewidth=0.5)
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Disorder threshold (0.5)')
    ax.axvline(x=aa_pos, color='darkred', linewidth=2, label=f'Truncation (aa {aa_pos})')
    ax.axvspan(aa_pos, len(sequence), alpha=0.15, color='red', label='Removed region')
    ax.set_ylabel('Disorder Score')
    ax.set_title(f'{gene} \u2014 Intrinsic Disorder Profile', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(0, 1)

    hydro = np.array([HYDROPHOBICITY.get(aa, 0) for aa in sequence])
    hydro_smooth = np.convolve(hydro, np.ones(window)/window, mode='same')
    ax = axes[1]
    ax.fill_between(positions, hydro_smooth, where=hydro_smooth >= 0,
                    alpha=0.3, color='orange', label='Hydrophobic')
    ax.fill_between(positions, hydro_smooth, where=hydro_smooth < 0,
                    alpha=0.3, color='blue', label='Hydrophilic')
    ax.plot(positions, hydro_smooth, color='gray', linewidth=0.5)
    ax.axvline(x=aa_pos, color='darkred', linewidth=2)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel('Hydrophobicity\n(K-D)')
    ax.legend(loc='upper right', fontsize=8)

    charges = np.array([CHARGE.get(aa, 0) for aa in sequence])
    charge_smooth = np.convolve(charges, np.ones(window)/window, mode='same')
    ax = axes[2]
    ax.fill_between(positions, charge_smooth, where=charge_smooth >= 0,
                    alpha=0.3, color='blue', label='Positive')
    ax.fill_between(positions, charge_smooth, where=charge_smooth < 0,
                    alpha=0.3, color='red', label='Negative')
    ax.plot(positions, charge_smooth, color='gray', linewidth=0.5)
    ax.axvline(x=aa_pos, color='darkred', linewidth=2)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel('Charge Density')
    ax.set_xlabel('Amino Acid Position')
    ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    outpath = os.path.join(output_dir, f'{gene}_disorder_profile.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved plot: {outpath}")


def plot_summary_comparison(results_df, output_dir):
    if len(results_df) == 0:
        return
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    genes = results_df['gene'].values
    x = np.arange(len(genes))
    for ax, col_rm, col_ret, ylabel, title in [
        (axes[0, 0], 'removed_mean_disorder', 'retained_mean_disorder',
         'Mean Disorder Score', 'Disorder: Removed vs Retained'),
        (axes[0, 1], 'removed_charge_per_residue', 'retained_charge_per_residue',
         'Charge per Residue', 'Charge: Removed vs Retained'),
        (axes[1, 0], 'removed_frac_disordered', 'retained_frac_disordered',
         'Fraction Disordered (>0.5)', 'Fraction Disordered: Removed vs Retained'),
        (axes[1, 1], 'removed_mean_hydrophobicity', 'retained_mean_hydrophobicity',
         'Mean Hydrophobicity (K-D)', 'Hydrophobicity: Removed vs Retained'),
    ]:
        ax.bar(x - 0.2, results_df[col_rm], 0.35, label='Removed region',
               color='indianred', alpha=0.8)
        ax.bar(x + 0.2, results_df[col_ret], 0.35, label='Retained protein',
               color='steelblue', alpha=0.8)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
    plt.suptitle('Truncation Analysis: Properties of Removed vs Retained Regions',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    outpath = os.path.join(output_dir, 'truncation_summary_comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved summary plot: {outpath}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Disorder analysis of truncation regions in NMD-escaping variants')
    parser.add_argument('--hc-variants', required=True,
                        help='Path to high_confidence_variants.xlsx')
    parser.add_argument('--output-dir', default='./disorder_results/',
                        help='Output directory (default: ./disorder_results/)')
    parser.add_argument('--min-position-fraction', type=float, default=0.0,
                        help='Only analyze variants with position_fraction >= this '
                             '(e.g., 0.85 for near-C-terminal only). Default: 0.0 (all)')
    parser.add_argument('--fetch-fasta', metavar='OUTPUT_FASTA', default=None,
                        help='If set, fetch UniProt sequences for all variants and '
                             'write to this FASTA file, then exit. No disorder '
                             'analysis is run. Output is ready for make_parsev2_fasta.py.')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --fetch-fasta mode: just get sequences and exit
    if args.fetch_fasta:
        fetch_fasta_only(args.hc_variants, args.fetch_fasta, args.min_position_fraction)
        sys.exit(0)

    # --- full disorder analysis below (unchanged) ---

    print("=" * 70)
    print("Intrinsic Disorder Analysis of NMD-Escaping Truncation Variants")
    print("=" * 70)

    print(f"\nStep 1: Loading variants from {args.hc_variants}...")
    variants = load_variants(args.hc_variants, args.min_position_fraction)
    if not variants:
        print("ERROR: No variants loaded. Check file path and filters.")
        sys.exit(1)
    for v in variants:
        print(f"  {v['gene']:12s} aa {v['aa_pos']}/{v['prot_len']} "
              f"({v['position_fraction']:.1%})  [{v['variant_type']}]")

    print(f"\nStep 2: Looking up UniProt accessions and fetching sequences...")
    gene_cache = {}
    all_results = []

    for variant in variants:
        gene = variant['gene']
        aa_pos = variant['aa_pos']
        prot_len = variant['prot_len']

        print(f"\n--- {gene} ({variant['variant_id']}) ---")
        print(f"  Truncation: aa {aa_pos} / {prot_len} "
              f"({variant['position_fraction']:.1%}), "
              f"removing {prot_len - aa_pos} residues")

        if gene in gene_cache:
            uniprot_id, sequence, header = gene_cache[gene]
            print(f"  UniProt: {uniprot_id} (cached)")
        else:
            print(f"  Looking up UniProt accession for {gene}...")
            uniprot_id, uniprot_len = lookup_uniprot_id(gene)
            if uniprot_id is None:
                print(f"  SKIPPING: Could not find UniProt accession")
                continue
            print(f"  Found: {uniprot_id} (length: {uniprot_len})")
            time.sleep(0.5)
            print(f"  Fetching sequence...")
            sequence, header = fetch_uniprot_sequence(uniprot_id)
            if sequence is None:
                print(f"  SKIPPING: Could not fetch sequence")
                continue
            gene_cache[gene] = (uniprot_id, sequence, header)
            time.sleep(0.5)

        actual_len = len(sequence)
        if actual_len != prot_len:
            print(f"  NOTE: UniProt length ({actual_len}) vs expected ({prot_len})")

        if aa_pos >= actual_len:
            print(f"  WARNING: Truncation position ({aa_pos}) >= sequence length "
                  f"({actual_len}), skipping")
            continue

        print(f"  Predicting disorder ({actual_len} residues)...")
        disorder_scores = predict_disorder(sequence)

        retained_seq = sequence[:aa_pos]
        removed_seq = sequence[aa_pos:]
        retained_disorder = disorder_scores[:aa_pos]
        removed_disorder = disorder_scores[aa_pos:]

        retained_props = analyze_region(retained_seq)
        removed_props = analyze_region(removed_seq)

        retained_mean_dis = np.mean(retained_disorder) if len(retained_disorder) > 0 else np.nan
        removed_mean_dis = np.mean(removed_disorder) if len(removed_disorder) > 0 else np.nan
        retained_frac_dis = np.mean(retained_disorder > 0.5) if len(retained_disorder) > 0 else np.nan
        removed_frac_dis = np.mean(removed_disorder > 0.5) if len(removed_disorder) > 0 else np.nan

        tail_50 = disorder_scores[-50:] if actual_len >= 50 else disorder_scores
        tail_50_mean = np.mean(tail_50)
        tail_50_frac_dis = np.mean(tail_50 > 0.5)

        flags = []
        if removed_mean_dis > 0.5 or removed_frac_dis > 0.5:
            flags.append("REMOVED_REGION_IS_IDR")
        if removed_props['net_charge'] < -1:
            flags.append("REMOVED_REGION_IS_ACIDIC")
        if removed_props['frac_proline'] + removed_props['frac_glycine'] > 0.2:
            flags.append("HIGH_PRO_GLY")
        if removed_props['frac_aromatic'] > 0.05:
            flags.append("AROMATIC_STICKERS_LOST")

        result = {
            'gene': gene,
            'variant_id': variant['variant_id'],
            'variant_type': variant['variant_type'],
            'uniprot_accession': uniprot_id,
            'uniprot_length': actual_len,
            'truncation_position': aa_pos,
            'protein_length_expected': prot_len,
            'removed_residues': actual_len - aa_pos,
            'position_fraction': aa_pos / actual_len,
            'pooled_weighted_efficiency': variant['pooled_weighted_efficiency'],
            'tier': variant['tier'],
            'gnomad_pLI': variant['gnomad_pLI'],
            'gnomad_LOEUF': variant['gnomad_LOEUF'],
            'retained_mean_disorder': retained_mean_dis,
            'removed_mean_disorder': removed_mean_dis,
            'retained_frac_disordered': retained_frac_dis,
            'removed_frac_disordered': removed_frac_dis,
            'tail_50_mean_disorder': tail_50_mean,
            'tail_50_frac_disordered': tail_50_frac_dis,
            'retained_mean_hydrophobicity': retained_props['mean_hydrophobicity'],
            'removed_mean_hydrophobicity': removed_props['mean_hydrophobicity'],
            'retained_net_charge': retained_props['net_charge'],
            'removed_net_charge': removed_props['net_charge'],
            'retained_charge_per_residue': retained_props['charge_per_residue'],
            'removed_charge_per_residue': removed_props['charge_per_residue'],
            'removed_frac_charged': removed_props['frac_charged'],
            'removed_frac_positive': removed_props['frac_positive'],
            'removed_frac_negative': removed_props['frac_negative'],
            'removed_frac_aromatic': removed_props['frac_aromatic'],
            'removed_frac_proline': removed_props['frac_proline'],
            'removed_frac_glycine': removed_props['frac_glycine'],
            'condensate_flags': '; '.join(flags) if flags else 'none',
            'removed_sequence': removed_seq,
        }
        all_results.append(result)

        print(f"  Retained (1-{aa_pos}):  disorder={retained_mean_dis:.3f}  "
              f"frac_dis={retained_frac_dis:.3f}  "
              f"charge/res={retained_props['charge_per_residue']:.4f}")
        print(f"  Removed  ({aa_pos+1}-{actual_len}): disorder={removed_mean_dis:.3f}  "
              f"frac_dis={removed_frac_dis:.3f}  "
              f"charge/res={removed_props['charge_per_residue']:.4f}  "
              f"net_charge={removed_props['net_charge']:.1f}")
        print(f"  Condensate flags: {result['condensate_flags']}")

        plot_disorder_profile(gene, sequence, disorder_scores, aa_pos, args.output_dir)

    if not all_results:
        print("\nERROR: No results generated. Check UniProt connectivity.")
        sys.exit(1)

    results_df = pd.DataFrame(all_results)
    csv_cols = [c for c in results_df.columns if c != 'removed_sequence']
    results_file = os.path.join(args.output_dir, 'disorder_analysis_results.csv')
    results_df[csv_cols].to_csv(results_file, index=False)
    print(f"\nSaved: {results_file}")

    seq_file = os.path.join(args.output_dir, 'removed_sequences.fasta')
    with open(seq_file, 'w') as f:
        for _, row in results_df.iterrows():
            f.write(f">{row['gene']}|{row['uniprot_accession']}|"
                    f"removed_aa{row['truncation_position']+1}-{row['uniprot_length']}|"
                    f"{row['variant_id']}\n")
            f.write(f"{row['removed_sequence']}\n")
    print(f"Saved: {seq_file}")

    print("\nGenerating summary plots...")
    plot_summary_comparison(results_df, args.output_dir)

    print("\n" + "=" * 70)
    print("DISORDER ANALYSIS SUMMARY")
    print("=" * 70)
    print(f"\nVariants analyzed: {len(results_df)}")
    print(f"Disorder predictor: {'metapredict' if USE_METAPREDICT else 'amino acid propensity (fallback)'}")

    print(f"\n{'Gene':10s} {'UniProt':8s} {'Trunc':>6s} {'Remvd':>6s} "
          f"{'Dis_rm':>8s} {'Dis_ret':>8s} {'Chg_rm':>8s} {'NMD_eff':>8s} {'Flags'}")
    print("-" * 100)
    for _, row in results_df.iterrows():
        print(f"{row['gene']:10s} {row['uniprot_accession']:8s} "
              f"{row['truncation_position']:6d} {row['removed_residues']:6d} "
              f"{row['removed_mean_disorder']:8.3f} {row['retained_mean_disorder']:8.3f} "
              f"{row['removed_charge_per_residue']:8.4f} "
              f"{row['pooled_weighted_efficiency']:8.3f} "
              f"{row['condensate_flags']}")

    n_idr = results_df['condensate_flags'].str.contains('REMOVED_REGION_IS_IDR').sum()
    n_acidic = results_df['condensate_flags'].str.contains('REMOVED_REGION_IS_ACIDIC').sum()
    print(f"\nCondensate-relevant findings:")
    print(f"  Variants where removed region is IDR:    {n_idr} / {len(results_df)}")
    print(f"  Variants where removed region is acidic: {n_acidic} / {len(results_df)}")

    print(f"\n{'='*70}")
    print(f"Output files in: {args.output_dir}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()