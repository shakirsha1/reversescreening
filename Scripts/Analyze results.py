#!/usr/bin/env python3
"""
Analyze PharmacoNet Reverse Screening Results
=============================================

Generate visualizations and statistical reports from screening results.

Author: PharmacoNet Team
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(
        description="Analyze reverse screening results"
    )
    parser.add_argument(
        "results_csv", type=str, help="screening_results.csv file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="analysis",
        help="output directory for plots and reports",
    )
    args = parser.parse_args()

    # Read results
    results_path = Path(args.results_csv)
    if not results_path.exists():
        print(f"ERROR: Results file not found: {args.results_csv}")
        sys.exit(1)

    print(f"Loading results from: {args.results_csv}")
    df = pd.read_csv(args.results_csv)

    if df.empty:
        print("WARNING: Results file is empty!")
        sys.exit(0)

    print(f"Loaded {len(df)} screening results")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate analysis report
    report_path = output_dir / "analysis_report.txt"
    print(f"Generating report: {report_path}")

    with open(report_path, "w") as f:
        f.write("="*80 + "\n")
        f.write("PharmacoNet Reverse Screening Analysis Report\n")
        f.write("="*80 + "\n\n")

        f.write("SUMMARY STATISTICS\n")
        f.write("-"*80 + "\n")
        f.write(f"Total matches: {len(df)}\n")
        f.write(f"Unique queries: {df['query_name'].nunique()}\n")
        f.write(f"Unique targets: {df['pharmacophore_model'].nunique()}\n")
        f.write(f"\nScore Statistics:\n")
        f.write(f"  Mean:   {df['score'].mean():.4f}\n")
        f.write(f"  Median: {df['score'].median():.4f}\n")
        f.write(f"  Std:    {df['score'].std():.4f}\n")
        f.write(f"  Min:    {df['score'].min():.4f}\n")
        f.write(f"  Max:    {df['score'].max():.4f}\n")

        f.write(f"\n\nTOP 10 QUERY-TARGET PAIRS\n")
        f.write("-"*80 + "\n")
        top_10 = df.nlargest(10, 'score')
        for idx, row in top_10.iterrows():
            model_name = Path(row['pharmacophore_model']).stem
            f.write(f"{row['query_name']:20s} → {model_name:40s} Score: {row['score']:.4f}\n")

        f.write(f"\n\nPER-QUERY SUMMARY\n")
        f.write("-"*80 + "\n")
        query_stats = df.groupby('query_name')['score'].agg(['count', 'mean', 'max'])
        f.write(query_stats.to_string())

    print("Report saved ✓")

    # Generate visualizations
    sns.set_style("whitegrid")

    # 1. Score distribution histogram
    print("Generating score distribution plot...")
    plt.figure(figsize=(10, 6))
    plt.hist(df['score'], bins=50, edgecolor='black', alpha=0.7)
    plt.xlabel('Pharmacophore Score', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Score Distribution Across All Matches', fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "score_distribution.png", dpi=300)
    plt.close()
    print("  → score_distribution.png saved ✓")

    # 2. Per-query boxplot
    if df['query_name'].nunique() > 1:
        print("Generating per-query boxplot...")
        plt.figure(figsize=(12, 6))
        df_sorted = df.sort_values('query_name')
        sns.boxplot(data=df_sorted, x='query_name', y='score')
        plt.xlabel('Query Molecule', fontsize=12)
        plt.ylabel('Score', fontsize=12)
        plt.title('Score Distribution by Query Molecule', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / "score_by_query_boxplot.png", dpi=300)
        plt.close()
        print("  → score_by_query_boxplot.png saved ✓")

    # 3. Heatmap (if multiple queries and targets)
    if df['query_name'].nunique() > 1 and df['pharmacophore_model'].nunique() > 5:
        print("Generating heatmap...")
        
        # Get top targets per query for visualization
        top_targets = df.groupby('query_name').apply(
            lambda x: x.nlargest(10, 'score')
        ).reset_index(drop=True)
        
        # Create pivot table
        pivot = top_targets.pivot_table(
            index='query_name',
            columns='pharmacophore_model',
            values='score',
            fill_value=0
        )
        
        # Simplify model names
        pivot.columns = [Path(col).stem[:30] for col in pivot.columns]
        
        plt.figure(figsize=(16, 8))
        sns.heatmap(
            pivot,
            cmap='YlOrRd',
            cbar_kws={'label': 'Pharmacophore Score'},
            linewidths=0.5,
            linecolor='gray'
        )
        plt.xlabel('Target Protein', fontsize=12)
        plt.ylabel('Query Molecule', fontsize=12)
        plt.title('Query-Target Score Heatmap (Top 10 per query)', fontsize=14, fontweight='bold')
        plt.xticks(rotation=90, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=10)
        plt.tight_layout()
        plt.savefig(output_dir / "query_target_heatmap.png", dpi=300)
        plt.close()
        print("  → query_target_heatmap.png saved ✓")

    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE ✓")
    print(f"{'='*80}")
    print(f"Output directory: {output_dir}")
    print(f"Files generated:")
    print(f"  - analysis_report.txt")
    print(f"  - score_distribution.png")
    if df['query_name'].nunique() > 1:
        print(f"  - score_by_query_boxplot.png")
    if df['query_name'].nunique() > 1 and df['pharmacophore_model'].nunique() > 5:
        print(f"  - query_target_heatmap.png")


if __name__ == "__main__":
    main()
