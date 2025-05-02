import pandas as pd
import numpy as np

def compute_bin_baf(rdr_file, baf_file, output_file=None):
    """
    Compute bin-level BAF values from SNP-level data using optimized pandas operations.
    
    Args:
        rdr_file (str): Path to the RDR file
        baf_file (str): Path to the BAF file
        output_file (str, optional): Path to save output TSV file
        
    Returns:
        pandas.DataFrame: DataFrame with bin-level BAF values
    """
    try:
        # Load data with appropriate columns and data types
        rdr_columns = ['CHROMOSOME', 'START', 'END', 'CELL', 'NORMAL', 'COUNT', 'RDR']
        rdr_df = pd.read_csv(rdr_file, sep='\t', header=None, names=rdr_columns)
        rdr_df = rdr_df[['CHROMOSOME', 'START', 'END', 'CELL']]  # Keep only relevant columns
        
        baf_columns = ['CHROMOSOME', 'POS', 'CELL', 'A_COUNT', 'B_COUNT']
        baf_df = pd.read_csv(baf_file, sep='\t', header=None, names=baf_columns)
        
        # Check for empty data
        if rdr_df.empty or baf_df.empty:
            print("Warning: Empty input data")
            return pd.DataFrame(columns=['CHROMOSOME', 'START', 'END', 'CELL', 'BAF'])
        
        # Convert to appropriate dtypes
        rdr_df['START'] = rdr_df['START'].astype(int)
        rdr_df['END'] = rdr_df['END'].astype(int)
        baf_df['POS'] = baf_df['POS'].astype(int)
        
        # Filter out SNPs with zero total counts
        baf_df = baf_df[(baf_df['A_COUNT'] + baf_df['B_COUNT']) > 0].copy()
        if baf_df.empty:
            print("Warning: No SNPs with valid counts")
            return pd.DataFrame(columns=['CHROMOSOME', 'START', 'END', 'CELL', 'BAF'])
        
      # Calculate min counts and total counts
        baf_df['MIN_COUNT'] = np.minimum(baf_df['A_COUNT'], baf_df['B_COUNT'])
        baf_df['TOTAL_COUNT'] = baf_df['A_COUNT'] + baf_df['B_COUNT']
        
        # Initialize results dataframe
        result_rows = []
        
        # Process each chromosome separately
        for chrom in rdr_df['CHROMOSOME'].unique():
            chrom_bins = rdr_df[rdr_df['CHROMOSOME'] == chrom].copy()
            chrom_snps = baf_df[baf_df['CHROMOSOME'] == chrom].copy()
            
            if chrom_snps.empty:
                # No SNPs for this chromosome, add all bins with NaN BAFs
                for _, bin_row in chrom_bins.iterrows():
                    result_rows.append({
                        'CHROMOSOME': bin_row['CHROMOSOME'],
                        'START': bin_row['START'],
                        'END': bin_row['END'],
                        'CELL': bin_row['CELL'],
                        'BAF': np.nan
                    })
                continue
            
            # Process each cell separately within this chromosome
            for cell in chrom_bins['CELL'].unique():
                cell_bins = chrom_bins[chrom_bins['CELL'] == cell].copy()
                cell_snps = chrom_snps[chrom_snps['CELL'] == cell].copy()
                
                if cell_snps.empty:
                    # No SNPs for this cell, add all bins with NaN BAFs
                    for _, bin_row in cell_bins.iterrows():
                        result_rows.append({
                            'CHROMOSOME': bin_row['CHROMOSOME'],
                            'START': bin_row['START'],
                            'END': bin_row['END'],
                            'CELL': bin_row['CELL'],
                            'BAF': np.nan
                        })
                    continue
                
                # Create bin edges for pd.cut (must include the rightmost edge)
                bin_edges = sorted(list(cell_bins['START']) + [max(cell_bins['END'])])
                
                # Assign each SNP to its bin
                cell_snps.loc[:, 'BIN_IDX'] = pd.cut(
                    cell_snps['POS'], 
                    bins=bin_edges,
                    labels=False,
                    include_lowest=True
                )
                
                # Group SNPs by bin and calculate sums
                bin_stats = cell_snps.groupby('BIN_IDX').agg({
                    'MIN_COUNT': 'sum',
                    'TOTAL_COUNT': 'sum'
                }).reset_index()
                
                # Convert bin index to bin range for joining
                bin_stats['BIN_START'] = [bin_edges[int(idx)] for idx in bin_stats['BIN_IDX']]
                
                # Create a mapping from bin start to BAF
                baf_map = {}
                for _, row in bin_stats.iterrows():
                    bin_start = row['BIN_START']
                    min_sum = row['MIN_COUNT']
                    total_sum = row['TOTAL_COUNT']
                    baf_map[bin_start] = min_sum / total_sum if total_sum > 0 else np.nan
                
                # Add all bins for this cell with calculated BAFs
                for _, bin_row in cell_bins.iterrows():
                    bin_start = bin_row['START']
                    baf = baf_map.get(bin_start, np.nan)
                    
                    result_rows.append({
                        'CHROMOSOME': bin_row['CHROMOSOME'],
                        'START': bin_row['START'],
                        'END': bin_row['END'],
                        'CELL': bin_row['CELL'],
                        'BAF': baf
                    })
        
        # Create final dataframe
        final_df = pd.DataFrame(result_rows)
        
        # Save output if requested
        if output_file:
            final_df.to_csv(output_file, sep='\t', index=False, header=False)
            print(f"Results saved to {output_file}")
        
        return final_df
    
    except Exception as e:
        print(f"Error processing data: {str(e)}")
        return pd.DataFrame(columns=['CHROMOSOME', 'START', 'END', 'CELL', 'BAF'])


rdr_file = '/alina-data1/Ezekiel/NSF_Cancer/chisel/src/rdr/rdr.tsv'
baf_file = '/alina-data1/Ezekiel/NSF_Cancer/chisel/src/baf/baf.tsv'
output_file = '/alina-data1/Ezekiel/NSF_Cancer/bins.tsv'
bin_df = compute_bin_baf(rdr_file, baf_file)