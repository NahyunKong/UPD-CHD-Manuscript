#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description='Merge multiple AnnotSV TSVs and filter variants.')
    parser.add_argument(
        'input_files', nargs='+',
        help='Paths to input TSV files.')
    parser.add_argument(
        '-o', '--output', required=True,
        help='Path to the output filtered TSV file.')
    args = parser.parse_args()

    # Load and concatenate all input files
    dfs = []
    for f in args.input_files:
        df = pd.read_csv(f, sep='\t', dtype=str)
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)

    # Ensure numeric columns are correctly typed
    combined['Gene_count'] = pd.to_numeric(
        combined.get('Gene_count', 0), errors='coerce').fillna(0).astype(int)

    # Step 2: Remove variants where RE_gene is empty and Gene_count == 0
    mask_no_gene = (
        combined.get('RE_gene', '').isna() |
        (combined['RE_gene'].str.strip() == '')
    ) & (combined['Gene_count'] == 0)
    # Create a copy to avoid SettingWithCopyWarning
    filtered = combined.loc[~mask_no_gene].copy()

    # Step 3: Remove variants with low ACMG or FULL classifications
    # Convert classification columns to numeric if they exist, using .loc to avoid warnings
    if 'ACMG_class' in filtered.columns:
        filtered.loc[:, 'ACMG_class'] = pd.to_numeric(
            filtered['ACMG_class'], errors='coerce')
    if 'full' in filtered.columns:
        filtered.loc[:, 'full'] = pd.to_numeric(
            filtered['full'], errors='coerce')

    mask_acmg_full = pd.Series(False, index=filtered.index)
    if 'ACMG_class' in filtered.columns:
        mask_acmg_full |= filtered['ACMG_class'].isin([1, 2])
    if 'full' in filtered.columns:
        mask_acmg_full |= filtered['full'].isin([1, 2])
    filtered = filtered.loc[~mask_acmg_full]

    # Write out the filtered table
    filtered.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
